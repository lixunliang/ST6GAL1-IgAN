const fs=require("fs");
const {Document,Packer,Paragraph,TextRun,Table,TableRow,TableCell,WidthType,BorderStyle,AlignmentType,PageOrientation,Footer,PageNumber}=require("docx");
const T=JSON.parse(fs.readFileSync("supp/tables.json","utf8"));
const F="Times New Roman";
function runs(t,opt={}){const out=[];const re=/10\^(-?\d+)/g;let last=0,m;
  while((m=re.exec(t))){out.push(new TextRun({text:t.slice(last,m.index)+"10",font:F,...opt}));out.push(new TextRun({text:m[1].replace("-","−"),font:F,superScript:true,...opt}));last=re.lastIndex;}
  out.push(new TextRun({text:t.slice(last),font:F,...opt}));return out;}
const W=12960; // landscape letter text width (15840-2*1440)
const NB={style:BorderStyle.NONE,size:0,color:"FFFFFF"};
const sections=[];
const note=fs.readFileSync("supp/supp_note1.md","utf8").split(/\n/).filter(l=>l.trim());
const nk=note.map((l,i)=>{ if(l.startsWith("# ")) return new Paragraph({children:runs(l.slice(2),{size:24,bold:true}),spacing:{after:240,line:480}});
  const m=l.match(/^([A-Z][^.]{2,40}\.) (.*)$/); 
  return new Paragraph({spacing:{line:480},children: m&&i>1? [...runs(m[1]+" ",{size:24,bold:true}),...runs(m[2],{size:24})] : runs(l,{size:24})}); });
sections.push({properties:{page:{size:{width:12240,height:15840},margin:{top:1440,bottom:1440,left:1440,right:1440}}},
  footers:{default:new Footer({children:[new Paragraph({alignment:AlignmentType.CENTER,children:[new TextRun({children:[PageNumber.CURRENT],font:F,size:20})]})]})},children:nk});
for(const [title,[head,rows,note]] of Object.entries(T)){
  const n=head.length; const cw=Array(n).fill(Math.floor(W/n)); cw[n-1]=W-cw[0]*(n-1);
  const mk=(cells,i,isHead)=>new TableRow({tableHeader:isHead,children:cells.map((c,j)=>new TableCell({width:{size:cw[j],type:WidthType.DXA},
     borders:{left:NB,right:NB,top:(isHead?{style:BorderStyle.SINGLE,size:12,color:"000000"}:NB),bottom:(isHead?{style:BorderStyle.SINGLE,size:6,color:"000000"}:(i===rows.length-1?{style:BorderStyle.SINGLE,size:12,color:"000000"}:NB))},
     margins:{top:30,bottom:30,left:60,right:60},
     children:[new Paragraph({children:runs(String(c),{size:16,bold:isHead})})]}))});
  const table=new Table({width:{size:W,type:WidthType.DXA},columnWidths:cw,rows:[mk(head,-1,true),...rows.map((r,i)=>mk(r,i,false))]});
  sections.push({properties:{page:{size:{width:12240,height:15840,orientation:PageOrientation.LANDSCAPE},margin:{top:1440,bottom:1440,left:1440,right:1440}}},
    footers:{default:new Footer({children:[new Paragraph({alignment:AlignmentType.CENTER,children:[new TextRun({children:[PageNumber.CURRENT],font:F,size:20})]})]})},
    children:[new Paragraph({children:runs(title,{size:22,bold:true}),spacing:{after:160}}),table,new Paragraph({children:runs(note,{size:18}),spacing:{before:120}})]});
}
const doc=new Document({styles:{default:{document:{run:{font:F,size:20}}}},sections});
Packer.toBuffer(doc).then(b=>fs.writeFileSync(process.argv[2]||"supp/Supplementary_Information.docx",b));
