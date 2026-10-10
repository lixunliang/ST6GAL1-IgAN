const fs=require("fs");const {Document,Packer,Paragraph,TextRun}=require("docx");
const F="Times New Roman";const L=fs.readFileSync("cover_letter.md","utf8").split(/\n/);
const kids=[];for(const l of L){kids.push(new Paragraph({spacing:{after:l.trim()?120:0,line:300},children:[new TextRun({text:l,font:F,size:24})]}));}
const doc=new Document({sections:[{properties:{page:{size:{width:12240,height:15840},margin:{top:1440,bottom:1440,left:1440,right:1440}}},children:kids}]});
Packer.toBuffer(doc).then(b=>fs.writeFileSync(process.argv[2],b));
