const fs=require("fs");
const {Document,Packer,Paragraph,TextRun,ImageRun,AlignmentType,Footer,PageNumber}=require("docx");
const md=fs.readFileSync("manuscript_v4.md","utf8");
const leg=md.split("## Figure legends")[1].split("## Supplementary")[0].split("\n").filter(l=>l.startsWith("Fig. "));
const imgs=["Fig1_design","Fig2_candidates","Fig3_two_signals","Fig4_celltype_phewas"];
const F="Times New Roman"; const kids=[];
imgs.forEach((f,i)=>{
  const buf=fs.readFileSync(f+".png"); const w=buf.readUInt32BE(16), h=buf.readUInt32BE(20);
  const W=620, H=Math.round(W*h/w);
  const [head,...rest]=leg[i].split(" | ");
  kids.push(new Paragraph({pageBreakBefore:i>0,alignment:AlignmentType.CENTER,children:[new ImageRun({type:"png",data:buf,transformation:{width:W,height:H}})]}));
  kids.push(new Paragraph({spacing:{before:240,line:480},children:[new TextRun({text:head+" | ",bold:true,font:F,size:24}),new TextRun({text:rest.join(" | "),font:F,size:24})]}));
});
const doc=new Document({styles:{default:{document:{run:{font:F,size:24}}}},sections:[{properties:{page:{size:{width:12240,height:15840},margin:{top:1440,bottom:1440,left:1440,right:1440}}},
 footers:{default:new Footer({children:[new Paragraph({alignment:AlignmentType.CENTER,children:[new TextRun({children:[PageNumber.CURRENT],font:F,size:24})]})]})},children:kids}]});
Packer.toBuffer(doc).then(b=>fs.writeFileSync(process.argv[2],b));
