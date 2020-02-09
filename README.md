# Visvariant

Abstract: Extremely large datasets are impossible or very difficult for humans to comprehend by standard mental approaches. Intuitive visualization of genetic variants in genomic sequencing data could help in the review and confirmation process of variants called by automated variant calling programs. To help facilitate interpretation of genetic variant next-generation sequencing (NGS) data we developed VisVariant, a customizable visualization tool that creates a figure showing the overlapping sequence information of thousands of individual reads including the variant and flanking regions.

# Example Usage:

java -Dfile.encoding=UTF-8 -Xmx4096m -classpath $Path/gv15OpenSource.jar: -jar $Path/gv15OpenSource.jar -datapath=$Path/ -referencepath=$Path/hg19_chr.fasta -variantpath=$Path/HG001_GRCh37_GIAB_chr18.vcf -phenotypepath=$Path/GIAB_phenotype2.csv -width=12000 -height=8000 -columns=30 -outputpath=$Path/ -r=chr18\;44595809\;A\;T

Where $Path is the linux full path,

Options:

-datapath=$Path/ # input bam files directory
-referencepath=$Path/hg19_chr.fasta # indexed reference file
-variantpath=$Path/HG001_GRCh37_GIAB_chr18.vcf # variants file
-phenotypepath=$Path/GIAB_phenotype2.csv # phenotype file
-outputpath=$Path/ # output directory
-r=chr18\;44595809\;A\;T # single snp 

# Run the program in parallel

my_func() {
cmd="java -Dfile.encoding=UTF-8 -Xmx4096m -classpath $Path/gv15OpenSource.jar: -jar $Path/gv15OpenSource.jar -datapath=$Path/ -referencepath=$Path/hg19_chr.fasta -variantpath=$Path/HG001_GRCh37_GIAB_chr18.vcf -phenotypepath=$Path/GIAB_phenotype2.csv -width=12000 -height=8000 -columns=30 -outputpath=$Path/ -r=$1"
echo $cmd
eval $cmd
}
export -f my_func

sed 's/\t/\\\;/g' HG001_GRCh37_GIAB_chr18.txt | parallel -j 4 my_func :::: -

# Common error where the program failed:

For example: SNP chr18\;6896450\;C\;T is failed, with the following output

Starting Engine
Processing variant 6896450
PacBio NoEdges = 60
Feb 09, 2020 8:03:00 AM gv15.Engine <init>
SEVERE: null
java.lang.NullPointerException
        at gv15.FragmentManager.ProcessFragments(FragmentManager.java:160)
        at gv15.Engine.<init>(Engine.java:132)
        at gv15.Gv15.main(Gv15.java:20)
        at sun.reflect.NativeMethodAccessorImpl.invoke0(Native Method)
        at sun.reflect.NativeMethodAccessorImpl.invoke(NativeMethodAccessorImpl.java:62)
        at sun.reflect.DelegatingMethodAccessorImpl.invoke(DelegatingMethodAccessorImpl.java:43)
        at java.lang.reflect.Method.invoke(Method.java:498)
        at com.sun.javafx.application.LauncherImpl.launchApplicationWithArgs(LauncherImpl.java:389)
        at com.sun.javafx.application.LauncherImpl.launchApplication(LauncherImpl.java:328)
        at sun.reflect.NativeMethodAccessorImpl.invoke0(Native Method)
        at sun.reflect.NativeMethodAccessorImpl.invoke(NativeMethodAccessorImpl.java:62)
        at sun.reflect.DelegatingMethodAccessorImpl.invoke(DelegatingMethodAccessorImpl.java:43)
        at java.lang.reflect.Method.invoke(Method.java:498)
        at sun.launcher.LauncherHelper$FXHelper.main(LauncherHelper.java:767)

Exam the bam file shows that not the whole region is covered by at least one reads
		
./samtools mpileup -f hg19_chr.fasta -r chr18:6896440-6896460 IonXpress_020_rawlib_chr18.b37.bam
[mpileup] 1 samples in 1 input files
chr18   6896446 G       65      ^u.^}.^}.^}.^}.^~.^}.^%.^~.^v.^%.^n.^%.^y.^~.^w.^%.^v.^x.^}.^{.^u.^~.^}.^y.^}.^p.^u.^}.^w.^}.^h.^w.^%.^~,^~,^~,^w,^~,^~,^~,^~,^~,^~,^~,^~,^~,^~,^~,^~,^~,^|,^~,^~,^~,^~,^~,^~,^~,^~,^~,^~,^~,^~,^~,  :<98=<<:=>;=7;?:;:;:4<5<:=9:<<;8;><<:8=;?><=B=>?><==?=?<>4<=<=9==
chr18   6896447 A       65      ..................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       <=<;==9==<==8<==8::<5=5<<?8;==<85==?/5>;9=<8<<<<A<><<<?=94<8<=/<<
chr18   6896448 G       65      ................................+1C..,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,    <<<;?=9=?<<;:;=<8::<8=:<;?8=<>;55==</5=:9=<8>=<;=;;=<:==96>8<=/9<
chr18   6896449 C       63      .................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, =<<;?=9@=<<;<==:8;:<9===;<8<<><5B9;/=4897889<7=5:<859;8485<=029
chr18   6896450 C       64      TTTTT.T..TTTTT.T.T.T.T.TTTTTT.T..,t,,t,t,,,,,t,t,tt,,,tt,t,ttt,,        <=<;B7969<::==9:6;6=4=5>;<?;;6<68=;98<9<>8====;C;:>>:<;</8@=<<9>
chr18   6896451 G       65      ..................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       <<=<><99><:<<===:<8?8<:=9<59;;:3;=<;95=9<=8=<==;C;<=>;=;=/8@=<>98
chr18   6896452 A       65      ..................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       <@==<<99=<<==<<=:<8=8=9?;>8;<=<5;>;>:5=7=>8?><<=@?<==<===/<=:<>=8
chr18   6896453 T       65      ..................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       >===<?9>>><<=<<=<=8=8=9=<>5<<;=8;=<<<5=8<><=G<<=A=<==<=<</<=:===8
chr18   6896454 A       65      ..................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       >===;>;=>>;<=;8=><<<;=;>>=:;<==5;8<<<:?8=<<=A=>=<=8?=<<==8<=;<=><
chr18   6896455 T       65      ..................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       >=<?;=;<?@;>=A8>@<=<==<=>=:<<=>5;8=895;59;79>9999929;:@9849998;98
chr18   6896456 T       65      ..................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       9989698899599929999969699968989/55@=>/?;??88>===@=8==<C=</>=><A<<
chr18   6896457 C       65      ..................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       ==<>;<<<>=;==>8<=<==;>;===:<><@5;9>==.>8>=68>==;=>6==<>=<7==>>=<<
chr18   6896458 A       65      ..................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       =<<=A<<<==<>>==<===>;=;=>=<=>>?5=>===/=;>>68?==<==6==<=<<7==C==<=
chr18   6896459 T       65      ................................+1C..,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,    >;<<==>8==9;=<<=B<=?==:<<>;<<8=5<C:=;7=;==<=>==>>=7=><><<:==<<<=>
chr18   6896460 C       64      .................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,        =;<<<?>8=>:;=<<=><>===8==>=<;8><>6975979999;99998/99797879999899


#References
  Multiple sequence alignments (MSA) visualization method - Sequence Bundles 
    http://bmcproc.biomedcentral.com/articles/10.1186/1753-6561-8-S2-S8 
    http://www.science-practice.com/blog/2015/08/25/sequence-bundles-web/
