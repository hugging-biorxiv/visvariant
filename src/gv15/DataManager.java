package gv15;

import com.opencsv.CSVReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ranasi01
 */
public class DataManager {
    
    private String filePath;
    private String variantPath;
    private String phenotypePath;
    private int phenotypeColumn;

    public DataManager(String filePath,String variantPath,String phenotypePath,int phenotypeColumn){
        this.filePath = filePath;
        this.variantPath = variantPath;
        this.phenotypePath = phenotypePath;
        this.phenotypeColumn = phenotypeColumn;
    }
    
    public File ImportVCFFile(){
        File vcfFile = new File(variantPath);
        if (!vcfFile.exists()) {
            System.out.println("ERROR: The VCF was not found. Check file location and address in prefs.txt file");
            System.exit(0);
        }
        return vcfFile;
    }

    public void ImportPhenotypes(HashMap<String,ArrayList<Phenotype>> phenotypes) {
        CSVReader reader = null;
        try {
            reader = new CSVReader(new FileReader(phenotypePath));
        } catch (FileNotFoundException ex) {
            System.out.println("ERROR: The CSV file was not found. Check file location and address in prefs.txt file");
            Logger.getLogger(DataManager.class.getName()).log(Level.SEVERE, null, ex);
            System.exit(0);
        }
        String[] nextLine;
        try {
            while ((nextLine = reader.readNext()) != null) { //@E look for the lines of the file
                // nextLine[] is an array of values from the line
                Phenotype tempPhenotype = new Phenotype();
                tempPhenotype.FID = nextLine[0]; //@E first term in the csv file
                tempPhenotype.IID = nextLine[1];
                tempPhenotype.FileName = nextLine[2]; //@E the DA00 part (in the test_regions)  
                
                if(!phenotypes.containsKey(nextLine[phenotypeColumn]))
                    phenotypes.put(nextLine[phenotypeColumn], new ArrayList());

                phenotypes.get(nextLine[phenotypeColumn]).add(tempPhenotype); 
            }
        } catch (IOException ex) {
            Logger.getLogger(DataManager.class.getName()).log(Level.SEVERE, null, ex);
        }
        phenotypes.remove("phenotype");
    }
    
    public void ImportPhenotypesOrder(ArrayList<String> phenotypesOrder) {
        CSVReader reader = null;
        boolean ElementExit = false;
        try {
            reader = new CSVReader(new FileReader(phenotypePath));
        } catch (FileNotFoundException ex) {
            System.out.println("ERROR: The CSV file was not found. Check file location and address in prefs.txt file");
            Logger.getLogger(DataManager.class.getName()).log(Level.SEVERE, null, ex);
            System.exit(0);
        }
        String[] nextLine;
        try {
            while ((nextLine = reader.readNext()) != null) {
                if (phenotypesOrder.isEmpty()) {
                    if (!nextLine[phenotypeColumn].equals("phenotype")) {
                        phenotypesOrder.add(nextLine[phenotypeColumn]);
                    }
                } else {
                    ElementExit = false;
                    for (int i = 0; i < phenotypesOrder.size(); i++) {
                        if (phenotypesOrder.get(i).equals(nextLine[phenotypeColumn])) {
                            ElementExit = true;
                        }
                    }
                    if (!ElementExit) {
                        phenotypesOrder.add(nextLine[phenotypeColumn]);
                    }
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(DataManager.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
