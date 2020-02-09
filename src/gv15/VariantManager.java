package gv15;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
/**
 *
 * @author ranasi01
 */
public class VariantManager {
    private Feature selectedVariant;
    private ArrayList<Feature> variants;
    final List<Allele> alleleList = new ArrayList<>();
    String lastName = "Lowe";
    public int TotalVariantCount = 0;
    
    public VariantManager(File vcfFile){
        
        variants = new ArrayList();        
        String variantCoordinate = UtilityFunctions.getInstance().VariantCoordinate;
        //If a specfic variant was provided, use the given coordinate ingore vcf file
        if(variantCoordinate!=null){
            String[] split = variantCoordinate.split(";");
            //  add ref allele
            alleleList.add(Allele.create(split[2], true));
            
            //  add one or more alt allele
            String[] alts = split[3].split(",");
            for (String alt : alts) {
                alleleList.add(Allele.create(alt, false));
            }
            
            // use Variant Context Builder to build a Variant Context variable
            VariantContextBuilder builder = new VariantContextBuilder("VQSR",split[0],
                    Integer.parseInt(split[1]),Integer.parseInt(split[1])+split[2].length() - 1,
                    alleleList);
            VariantContext vc =	builder.make();
            variants.add(vc);
            TotalVariantCount=1;
        }else{
            VCFFileReader reader = new VCFFileReader(vcfFile, false);
            for(VariantContext context:reader.iterator().toList()){ 
                variants.add(context);
                TotalVariantCount++;
            }            
        }

    }
    
    public Feature getSelectedVariant(){
        return selectedVariant;
    }
    
    public void setVariant(int index){
        selectedVariant = variants.get(index);
    }

    private void split(String string) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
