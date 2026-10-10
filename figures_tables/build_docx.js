const fs=require("fs");
const {Document,Packer,Paragraph,TextRun,HeadingLevel,AlignmentType,Footer,PageNumber,LineNumberRestartFormat,PageBreak}=require("docx");
const md=fs.readFileSync("manuscript_v4.md","utf8").split(/\n/);
const F="Times New Roman";
function runs(t,opt={}){
  const out=[]; const re=/10−(\d+)/g; let last=0,m;
  while((m=re.exec(t))){ out.push(new TextRun({text:t.slice(last,m.index)+"10",font:F,size:24,...opt})); out.push(new TextRun({text:"−"+m[1],font:F,size:24,superScript:true,...opt})); last=re.lastIndex; }
  out.push(new TextRun({text:t.slice(last),font:F,size:24,...opt})); return out;
}
const kids=[]; const sp={line:480,after:0};
const breakBefore=new Set(["Abstract","Introduction","References","Figure legends","Supplementary information (planned)"]);
for(const raw of md){
  const l=raw.trimEnd(); if(!l) continue;
  if(l.startsWith("# ")) kids.push(new Paragraph({children:runs(l.slice(2),{bold:true,size:28}),alignment:AlignmentType.LEFT,spacing:{line:480,after:240}}));
  else if(l.startsWith("## ")){ const h=l.slice(3); kids.push(new Paragraph({heading:HeadingLevel.HEADING_1,pageBreakBefore:breakBefore.has(h),children:runs(h,{bold:true,color:"000000"}),spacing:{line:480,before:240,after:0}})); }
  else if(l.startsWith("### ")) kids.push(new Paragraph({heading:HeadingLevel.HEADING_2,children:runs(l.slice(4),{bold:true,color:"000000"}),spacing:{line:480,before:120,after:0}}));
  else kids.push(new Paragraph({children:runs(l),spacing:sp}));
}
const doc=new Document({styles:{default:{document:{run:{font:F,size:24}}}},
 sections:[{properties:{page:{size:{width:12240,height:15840},margin:{top:1440,bottom:1440,left:1440,right:1440}},lineNumbers:{countBy:1,restart:LineNumberRestartFormat.CONTINUOUS}},
  footers:{default:new Footer({children:[new Paragraph({alignment:AlignmentType.CENTER,children:[new TextRun({children:[PageNumber.CURRENT],font:F,size:24})]})]})},
  children:kids}]});
Packer.toBuffer(doc).then(b=>fs.writeFileSync(process.argv[2]||"ST6GAL1_IgAN_manuscript_v3.docx",b));
