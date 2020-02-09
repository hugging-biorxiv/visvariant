package gv15;

/**
 *
 * @author ranasi01
 */
public class AppOptions {
    public String OptionsList = "Usage: gv15.jar [-options]\n"
            + "Available Options\n"
            + " -pref={FullPath}                Define the full path and file name of a parameters file. Other options are ignored when this is in use.\n"
            + " -reference={FullPath}           Define the full path and file name of the Reference .FASTA file.\n"
            + " -variant={FullPath}             Define the full path and file name of Variants .vcf file.\n"
            + " -phenotype={FullPath}           Define the full path and file name of the Phenotypes .csv file.\n"
            + " -outputpath={val}               Path to dump rendered images.\n"
            + " -datapath={val}                 Path of the target BAM and BAI files.\n"
            + " -r={val}                        Chromosome number and coordinate of target variant. Ignores vcf file if present. FORMAT -r=chr{chromosome_num}_{cooridnate}."
            + " -width={val}                    Width of the rendered image. DEFAULT -width=1900.\n"
            + " -height={val}                   Height of the rendered image. DEFAULT -height=960.\n"
            + " -flank={val}                    Number of base pairs to left/right of variant. DEFAULT -flank=15.\n"
            + " -gridstartx={val}               X coordinate of the first panel. DEFAULT -gridstartx=200.\n"
            + " -gridstarty={val}               Y coordinate of the first panel. DEFAULT -gridstarty=100.\n"
            + " -panelseparation={val}          Seperation between panels in pixels. DEFAULT -panelseparation=215.\n"
            + " -columns={val}                  Number of columns to render in each panel. DEFAULT -columns=37.\n"
            + " -columnwidth={val}              Width of each individual column in all panels. DEFAULT -columnwidth=45.\n"
            + " -rowheight={val}                Height of each individual row in all panels. DEFAULT -rowheight=18.\n"
            + " -fragmentxoffset={val}          Padding applied between fragments. DEFAULT -fragmentxoffset=15.\n"
            + " -outputtype={val}               Output type of rendered images. Supports png/jpeg. DEFAULT -outputtype=png.\n"
            + " -phenotypecolumn={val}          Target phenotype column in the phenotypes file. DEFAULT -phenotypecolumn=4.\n"
            + " -fontSizeMultiplier={val}       Enlarge fonts. DEFAULT -fontSizeMultiplier=4.\n"
            + " -lineThicknessMultiplier={val}  Enhance the lines. DEFAULT -lineThicknessMultiplier=4.\n"
            + " -readcountrenderthreshold={val} Minimum read count to render fragment. DEFAULT -readcountrenderthreshold=0.\n"
            + " -insertionsonlyatvariant={val}  Show all insertions 0. Show insertions only at variant 1. DEFAULT -insertionsonlyatvariant=0.\n"
            + " -readcolour_unvaried={val}      HEX colour of unvaried fragments. DEFAULT -readcolour_unvaried=#DCDCDC. Applicable colours https://docs.oracle.com/javafx/2/api/javafx/scene/paint/Color.html\n"
            + " -readcolour_varied={val}        HEX colour of varies fragments. DEFAULT -readcolour_varied=#FFA500. Applicable colours https://docs.oracle.com/javafx/2/api/javafx/scene/paint/Color.html\n"
            + " -readcolour_insertion={val}     HEX colour of insertion fragments. DEFAULT -readcolour_insertion=#8A2BE2. Applicable colours https://docs.oracle.com/javafx/2/api/javafx/scene/paint/Color.html\n";
}
