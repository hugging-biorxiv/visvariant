package gv15;

import htsjdk.samtools.BAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.tribble.Feature;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ranasi01
 * @author eisaa
 */
public class ReadManagerSAMTOOLS {

    private HashMap<String, ArrayList<Read>> readCollection;
    private HashMap<String, ArrayList<InsertFeature>> InsertFeatures;
    private String dataPath;

    public HashMap<String, ArrayList<String>> loadedReferences;
    //public TabletDataHandler tabletDataHandler;
    public HashMap<String, int[]> InsertionArrays;
    public String referencePath;

    public ReadManagerSAMTOOLS(String referencePath, String dataPath) {
        this.referencePath = referencePath;
        readCollection = new HashMap(); // bam file stuff here
        //tabletDataHandler = new TabletDataHandler(cachePath); // remove this !!!
        loadedReferences = new HashMap(); // this is where the reference file is read in. Near the bottom. see if read in earlier elsewhere to compare. IT IS. read in earler and then use that read, instread of reading it new each time.
        this.dataPath = dataPath;
        InsertFeatures = new HashMap();
        InsertionArrays = new HashMap();
    }

    public void LoadDataFromSamples(HashMap<String, ArrayList<Phenotype>> phenotypes, // can use the phenotype
            int startCoordinate, int endCoordinate, ReferenceManager referenceManager, // can use the reference manager
            Feature currentVariant) throws FileNotFoundException, UnsupportedEncodingException {

        ReferenceSequence DNAsequence;
        File referenceFile = new File(referencePath); // READ Reference // must chnage this bit
        //FastaSequenceIndexCreator.create(referenceFilePath, true);         
        IndexedFastaSequenceFile reference;
        reference = new IndexedFastaSequenceFile(referenceFile);
        DNAsequence = reference.getSubsequenceAt(currentVariant.getContig(), startCoordinate, endCoordinate); // get the subsequence

        for (String type : phenotypes.keySet()) { // 4 atm (4 phenotypes) //THIS SECTION TAKES TIME

            ArrayList<String> referenceSequenceBases = new ArrayList<String>();
            byte[] byteArray = DNAsequence.getBases();
            String refSequenceInTargetRegion = new String(DNAsequence.getBaseString());
            for (int i = 0; i < refSequenceInTargetRegion.length(); i++) {
                referenceSequenceBases.add(refSequenceInTargetRegion.charAt(i) + "");
            }

            loadedReferences.put(type, referenceSequenceBases); // what in in the loaded references
            referenceManager.AddReference(type, referenceSequenceBases);

            for (Phenotype currentSample : phenotypes.get(type)) { // each of the bam files

                //Sample and reference data input for Tablet
                String[] fileNames = new String[]{ // the bam and fasta file locations
                    dataPath + "/" + currentSample.FileName, referenceManager.getReferencePath()};

                SamReader samReader = SamReaderFactory.makeDefault(). // open the selected bam file
                        validationStringency(ValidationStringency.LENIENT).
                        open(new File(fileNames[0]));

                SAMRecordIterator BAMFileIterator = samReader.queryOverlapping(currentVariant.getContig().toString(),
                        startCoordinate, endCoordinate); // read all the bam files that cover the target region

                int count = 0; // diagnostics
                int sequenceLength;
                String sequenceRead;

                //Add the samples and its reads to the collection
                readCollection.put(currentSample.FileName, new ArrayList());
                while (BAMFileIterator.hasNext()) {
                    BAMRecord rec = (BAMRecord) BAMFileIterator.next();
                    count++;
                    sequenceLength = rec.getReadLength();

                    Cigar cigar = rec.getCigar();
                    Iterator<CigarElement> cigIter = cigar.iterator();
                    CigarElement currentCigElement = cigIter.next();
                    CigarOperator currentOperator = currentCigElement.getOperator();
                    int currentOperatorSize = currentCigElement.getLength();

                    ArrayList<String> cigarArray = new ArrayList<String>();
                    int portionOfCigarElementElapsed = 0;

                    int softClipCount = 0;
                    int insertionCount = 0;
                    int deletionCount = 0;

                    while (cigIter != null) {
                        if (portionOfCigarElementElapsed >= currentOperatorSize) {
                            if (cigIter.hasNext() == false) {
                                break;
                            }
                            currentCigElement = cigIter.next();
                            portionOfCigarElementElapsed = 0;
                            currentOperator = currentCigElement.getOperator(); // use the new cigar info
                            currentOperatorSize = currentCigElement.getLength();
                        }
                        //System.out.print(currentOperator);
                        cigarArray.add(currentOperator + "");
                        portionOfCigarElementElapsed++;
                    }

                    for (int i = 0; i < cigarArray.size(); i++) { // count the soft clippings and the insertions and deletions
                        if (cigarArray.get(i).equals("S")) {
                            softClipCount++;
                        }
                        if (cigarArray.get(i).equals("I")) {
                            insertionCount++;
                        }
                        if (cigarArray.get(i).equals("D")) {
                            deletionCount++;
                        }
                    }

                    sequenceRead = new String(rec.getReadBases(), "UTF-8");
                    Read tempRead = new Read();
                    tempRead.BaseValues = new String[sequenceRead.length() + deletionCount];
                    tempRead.StartPosition = rec.getAlignmentStart() - 1;


                    tempRead.Length = rec.getReadLength();
                    tempRead.ReadID = currentSample.toString() + "" + count;

                    int basePosition;
                    //for (int i = 0; i < sequenceLength - softClipCount - insertionCount + deletionCount; i++) { // the soft clips should nto be counted nor should the insertions
                    for (int i = 0; i <= rec.getAlignmentEnd() - rec.getAlignmentStart(); i++) { // the soft clips should nto be counted nor should the insertions    
                        
                        //basePosition = rec.getReadPositionAtReferencePosition(i + tempRead.StartPosition); // +1 because i starts from 0 but the read is 1 based
                        basePosition = getReadPositionFromReferenceZeroBased(i + tempRead.StartPosition, rec); // method to call the SAMtools fucntion dealing with 1-based data
                        if (basePosition == -1) { // there is a deletion
                            tempRead.BaseValues[i] = "N";
                        } else {
                            tempRead.BaseValues[i] = sequenceRead.charAt(basePosition) + "";
                        }
                    }

                    //Add the read to the collection
                    readCollection.get(currentSample.FileName).add(tempRead);

                    // inserts
                    int referencePosition;
                    int deletionOffset = 0; // deletions are in adjusted sequence but not original sequence so offset required
                    for (int i = 0; i < sequenceLength; i++) { // see what the base is at each postion and keep track of the actual location
                        if (cigarArray.get(i).equals("D")) {
                            deletionOffset++;
                            continue;
                        }

                        if (cigarArray.get(i).equals("I")) {
                            referencePosition = getReferencePositionFromReadZeroBased(i - deletionOffset - 1, rec); // get the postion the new bases are inserted at.
                            //-1 becuase we want to get the position of the base before the inserted base for to store the start position
                            
                            int currentInsertCount = getInsertSize(cigarArray, i);

                            if (referencePosition >= startCoordinate-1 && referencePosition < endCoordinate) { // if insert is in our area of interest
                                InsertFeature tempInsert = new InsertFeature();

                                tempInsert.TargetReadID = currentSample.toString() + "" + count;
                                tempInsert.StartCoodinate = referencePosition;

                                tempInsert.InsertedBases = new ArrayList<String>();

//                              System.out.println(currentInsertCount); // debug
                                tempInsert.InsertedBases = getInsertedBases(sequenceRead, i, currentInsertCount, deletionOffset);

                                if (!InsertFeatures.containsKey(currentSample.FileName)) {
                                    InsertFeatures.put(currentSample.FileName, new ArrayList()); //@E the first time add the file to the list
                                }
                                InsertFeatures.get(currentSample.FileName).add(tempInsert); //@E add the iserted base
                            }
                        }
                    }
                    count++; //thats one read
                }
                //System.out.println("file");
            }
        }
    }

    // These 2 methods are to make data 0-based
    public static int getReferencePositionFromReadZeroBased(int ReadPosition, BAMRecord record) {
        int referencePosition = record.getReferencePositionAtReadPosition(ReadPosition + 1); // to make the data 1-based for the input
        return referencePosition - 1; // make data 0-based
    }

    public static int getReadPositionFromReferenceZeroBased(int referencePosition, BAMRecord record) {
        int readPosition = record.getReadPositionAtReferencePosition(referencePosition + 1); // make the data 1-based
        return readPosition - 1; // make data 0 based            
    }

    public static int getInsertSize(ArrayList<String> cigarArray, int position) { // find number of inserts
        int count = 0;
        while (cigarArray.get(position).equals("I")) {
            count++;
            position++;
            if (position >= cigarArray.size()) {
                break;
            }
        }
        return count;
    }

    public static ArrayList<String> getInsertedBases(String sequenceReads, int position, int noOfInserts, int deletionOffset) {
        ArrayList<String> inserts = new ArrayList<String>();
        for (int i = 0; i < noOfInserts; i++) {
            inserts.add(sequenceReads.charAt(position - deletionOffset) + "");
            position++;
        }
        return inserts;
    }

    public void CreateInsertionArrays(HashMap<String, ArrayList<Phenotype>> phenotypes, int startCoordinate) {
        for (String type : phenotypes.keySet()) {
            //if(type.equals("Neg_Control")){

            InsertionArrays.put(type, new int[loadedReferences.get(type).size()]);

            for (Phenotype currentSample : phenotypes.get(type)) { // loop though all the files
                //If the current sample has any Inserts
                if (InsertFeatures.containsKey(currentSample.FileName)) {
                    for (InsertFeature feature : InsertFeatures.get(currentSample.FileName)) { // all the reads in a file // WE DO HAVE INSERT FEATURE

                        int insertedPos = (feature.StartCoodinate + 1) - startCoordinate;

                        if (insertedPos > -1) {

                            // if (feature.InsertedBases.size() == 6) { // why is this significant
                               // System.err.println("");
                            // }

                            if (InsertionArrays.get(type)[insertedPos] < feature.InsertedBases.size()) {
                                InsertionArrays.get(type)[insertedPos] = feature.InsertedBases.size();
                            }
                        }
                    }
                }
            }
            //}//End Type Check
        }
    }

    public ArrayList<Read> GetReadsForSample(String sampleName) { // do have read
        return readCollection.get(sampleName);
    }

    public int[] getInsertionArray(String phenotype) {
        return InsertionArrays.get(phenotype);
    }

    public ArrayList<InsertFeature> GetInsertsForReadAtPosition(Read read, String sampleName, // redo this
            int columnPosition, int startCoordinate) {
        ArrayList<InsertFeature> inserts = new ArrayList();
        if (InsertFeatures.containsKey(sampleName)) {
            for (InsertFeature feature : InsertFeatures.get(sampleName)) {
                if (feature.TargetReadID.equals(read.ReadID) && ((feature.StartCoodinate + 1) - startCoordinate) == columnPosition) {
                    inserts.add(feature);
                }
            }
        }

        return inserts;
    }

    private ArrayList<String> GetInsertedBases(String inserts) {
        ArrayList<String> bases = new ArrayList();

        for (int index = 0; index < inserts.length(); index++) {
            bases.add(Character.toString(inserts.charAt(index)));
        }

        return bases;
    }
}
