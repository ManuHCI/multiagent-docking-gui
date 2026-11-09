"""Local molecular views. No remote structures or CDN requests at runtime."""
from __future__ import annotations

import html
import json
import re
from pathlib import Path

from .project import Project


def values(project: Project, sid: str) -> dict:
    step = project.get_step(sid)
    if not step or step["status"] != "complete":
        return {}
    return json.loads(step["values_json"] or "{}")


def local_text(project: Project, filename: str) -> str:
    if not filename:
        return ""
    path = Path(filename).resolve()
    if not path.is_relative_to(project.folder.resolve()) or not path.is_file():
        return ""
    if path.stat().st_size > 8_000_000:
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def reference_atoms(pdb: str, label: str) -> str:
    parts = label.split()
    if len(parts) < 2:
        return ""
    match = re.fullmatch(r"(.)(-?\d+)([A-Za-z]?)", parts[1])
    if not match:
        return ""
    chain, number, insertion = match.groups()
    return "\n".join(line for line in pdb.splitlines()
                     if line.startswith("HETATM") and line[17:20].strip() == parts[0]
                     and line[21:22] == chain and line[22:26].strip() == number
                     and line[26:27].strip() == insertion
                     and line[76:78].strip().upper() not in {"H", "D"})


def payload(project: Project, stage_id: str, molecule: str = "", original: bool = False) -> dict:
    structure = values(project, "02_structure").get("chosen", {})
    receptor = values(project, "03_receptor")
    raw = local_text(project, structure.get("local_file", ""))
    prepared = local_text(project, receptor.get("prepared_file", ""))
    protein = raw if original or stage_id == "02_structure" else prepared or raw
    ref = reference_atoms(protein, receptor.get("ligand") or structure.get("best_ligand", ""))
    ref_lines = set(ref.splitlines())
    extras = "\n".join(line for line in protein.splitlines() if line.startswith("HETATM") and line not in ref_lines)
    # Keep all waters/cofactors out of the blue cartoon. The green reference
    # is loaded separately so its identity and colour are unambiguous.
    protein = "\n".join(line for line in protein.splitlines()
                        if line.startswith(("ATOM", "TER", "HELIX", "SHEET")))
    data = {"protein": protein, "reference": ref, "extras": extras, "ligand": "", "format": "pdb",
            "box": {}, "name": molecule, "focus": bool(ref), "overlay": False}
    if stage_id == "05_ligands":
        prepared_mols = values(project, "05_ligands").get("prepared", [])
        chosen = next((m for m in prepared_mols if m["name"] == molecule), next(iter(prepared_mols), {}))
        data.update(protein="", reference="", extras="", ligand=local_text(project, chosen.get("file", "")),
                    format="sdf", name=chosen.get("name", ""), focus=False)
    elif stage_id in {"04_site", "06_dock"}:
        data["box"] = values(project, "04_site").get("box", {})
        if stage_id == "06_dock":
            docking = values(project, "06_dock")
            results = docking.get("results", [])
            chosen = next((m for m in results if m.get("name") == molecule),
                          next((m for m in results if m.get("file")), {}))
            text = local_text(project, chosen.get("file", ""))
            if text:
                from stages.md import pdbqt_pose_to_pdb
                data["ligand"] = pdbqt_pose_to_pdb(text)
                data["name"] = chosen.get("name", "")
                # The recorded re-docking identity belongs to this target.
                # A familiar drug name alone never establishes reference identity.
                data["overlay"] = bool(ref and docking.get("redock_rmsd") is not None
                    and docking.get("redock_molecule")
                    and chosen.get("name", "").casefold() == docking["redock_molecule"].casefold())
    return data


def document(data: dict) -> str:
    encoded = json.dumps(data).replace("<", "\\u003c").replace("&", "\\u0026")
    return '''<!doctype html><html lang="en"><head><meta charset="utf-8">
<style>
body{margin:0;background:#f6f9fc;font:14px system-ui;color:#21354b;display:flex;flex-direction:column;height:100vh}
#view{width:100%;flex:1;min-height:200px;position:relative}
.tools{display:flex;gap:8px;padding:10px;flex-wrap:wrap;flex-shrink:0;border-top:1px solid #dce5ec}
button{background:white;border:1px solid #c3d3df;border-radius:7px;padding:8px 12px;color:#193d60;cursor:pointer}
button:focus-visible{outline:3px solid #1a73e8}#status{padding:8px 12px;font-size:12px;flex-shrink:0}
</style><script src="/assets/vendor/3Dmol-min.js"></script></head><body>
<div id="view" role="img" aria-label="Interactive molecular model. Drag to rotate and scroll to zoom."></div>
<div class="tools"><button id="whole">Show whole structure</button><button id="pocket">Focus on pocket</button>
<button id="box">Hide search box</button><button id="protein">Hide protein</button><button id="atoms">Show atoms</button></div>
<div id="status" role="status">Loading local 3D model…</div>
<script>
const data = ''' + encoded + ''';
const status = document.getElementById('status');
try {
 if (!window.$3Dmol) throw new Error('The bundled 3D library is missing.');
 const viewer = $3Dmol.createViewer(document.getElementById('view'), {backgroundColor:'#f6f9fc', antialias:true});
 let protein, reference, ligand, atomCount=0, atomsVisible=false;
 if (data.protein) { protein=viewer.addModel(data.protein,'pdb'); protein.setStyle({}, {cartoon:{color:'#6f9fcd',opacity:0.8}}); atomCount+=protein.selectedAtoms({}).length; }
 if (data.reference) { reference=viewer.addModel(data.reference,'pdb'); reference.setStyle({}, {stick:{color:'#219653',radius:0.2}}); atomCount+=reference.selectedAtoms({}).length; }
 if (data.extras) { const extras=viewer.addModel(data.extras,'pdb'); extras.setStyle({}, {stick:{colorscheme:'Jmol',radius:0.12}}); extras.setStyle({resn:['HOH','WAT','DOD']},{sphere:{radius:0.18,color:'#c65e65'}}); atomCount+=extras.selectedAtoms({}).length; }
 if (data.ligand) { ligand=viewer.addModel(data.ligand,data.format); ligand.setStyle({}, {stick:{colorscheme:'orangeCarbon',radius:0.2}}); atomCount+=ligand.selectedAtoms({}).length; }
 const b=data.box; let boxVisible=true, proteinVisible=true, lines=[];
 function drawBox() {
  lines.forEach(l=>viewer.removeShape(l)); lines=[];
  if(!boxVisible || !b || !b.size_x) return;
  const c=[b.center_x,b.center_y,b.center_z], s=[b.size_x,b.size_y,b.size_z];
  const pts=[]; for(let i=0;i<8;i++) pts.push({x:c[0]+((i&1)?1:-1)*s[0]/2,y:c[1]+((i&2)?1:-1)*s[1]/2,z:c[2]+((i&4)?1:-1)*s[2]/2});
  for(let i=0;i<8;i++) for(const bit of [1,2,4]) if(!(i&bit)) lines.push(viewer.addLine({start:pts[i],end:pts[i|bit],color:'#d28a00',linewidth:2}));
 }
 function focus() { if(reference) viewer.zoomTo({model:reference.getID()}); else if(ligand) viewer.zoomTo({model:ligand.getID()}); else viewer.zoomTo(); viewer.zoom(0.75); viewer.render(); }
 drawBox(); if(data.focus) focus(); else {viewer.zoomTo();viewer.render();}
 document.getElementById('whole').onclick=()=>{viewer.zoomTo();viewer.render();};
 document.getElementById('pocket').onclick=focus;
 document.getElementById('pocket').disabled=!reference&&!ligand;
 document.getElementById('box').hidden=!b||!b.size_x;
 document.getElementById('box').onclick=e=>{boxVisible=!boxVisible;drawBox();e.target.textContent=boxVisible?'Hide search box':'Show search box';viewer.render();};
 document.getElementById('protein').hidden=!protein;
 document.getElementById('protein').onclick=e=>{if(!protein)return;proteinVisible=!proteinVisible;proteinVisible?protein.show():protein.hide();e.target.textContent=proteinVisible?'Hide protein':'Show protein';viewer.render();};
 document.getElementById('atoms').onclick=e=>{if(!protein)return;atomsVisible=!atomsVisible;protein.setStyle({},atomsVisible?{stick:{colorscheme:'Jmol',radius:0.10}}:{cartoon:{color:'#6f9fcd',opacity:0.8}});e.target.textContent=atomsVisible?'Show ribbon':'Show atoms';viewer.render();};
 if(!protein) { document.getElementById('protein').hidden=true; document.getElementById('atoms').hidden=true; document.getElementById('pocket').textContent='Focus on molecule'; }
 new ResizeObserver(()=>{viewer.resize();viewer.render();}).observe(document.getElementById('view'));
 new IntersectionObserver(entries=>{if(entries.some(e=>e.isIntersecting)){viewer.resize();viewer.render();}}).observe(document.getElementById('view'));
 window.addEventListener('resize',()=>{viewer.resize();viewer.render();});
 if(!atomCount) throw new Error('No current structure is available. Run this step or load the saved example.');
 status.textContent='3D ready · '+atomCount+' atoms · Drag to rotate; scroll or pinch to zoom.';
} catch(error) { status.textContent='3D view unavailable: '+error.message+' Use the written results while checking WebGL support.'; }
</script></body></html>'''


def render(project: Project, stage_id: str, molecule: str = "", original: bool = False) -> None:
    from nicegui import ui
    data = payload(project, stage_id, molecule, original)
    if not data["protein"] and not data["ligand"]:
        ui.label("The 3D view appears when a current structure or prepared molecule is available.").classes("text-slate-500")
        return
    iframe = '<iframe title="Interactive molecular view" sandbox="allow-scripts" style="width:100%;height:445px;border:1px solid #dce5ec;border-radius:12px" srcdoc="' + html.escape(document(data), quote=True) + '"></iframe>'
    ui.html(iframe, sanitize=False).classes("w-full")
    ui.label("Blue: protein · Green: measured reference ligand · Orange: predicted/prepared molecule · Amber: search box").classes("text-xs text-slate-600")
    if stage_id == "06_dock" and data["ligand"]:
        ui.label(f"Reference overlay: verified same molecule ({data['name']})." if data["overlay"] else
                 "Same-molecule reference identity is not verified for this selection. Visual overlap is not a re-docking RMSD comparison.").classes("text-sm text-slate-600")
