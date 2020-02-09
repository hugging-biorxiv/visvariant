package gv15;


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.tribble.Feature;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author ranasi01
 * @author Eisa Anwar
 */
public class FragmentManager {
    //DebugInfo red = new DebugInfo();
    
    private int maxReadCount;
    private String dataPath;
    private ReadManagerSAMTOOLS readManager;
    
    public FragmentManager(String dataPath){
        this.dataPath = dataPath;
        this.maxReadCount = 0;
    }
    
    public void ProcessFragments(int j, HashMap<String,ArrayList<Phenotype>> phenotypes, ArrayList<String> phenotypesOrder, ReferenceManager referenceManager,
            PanelManager panelManager,int flank,Feature currentVariant) throws Exception{
        //@E j added sent from i
        int startCoord = currentVariant.getStart() - flank;
        int endCoord  = currentVariant.getStart() + flank;
        int loopCount =0;
        //Load the Read data for all Samples
        readManager = new ReadManagerSAMTOOLS(referenceManager.getReferencePath(),dataPath);
        readManager.LoadDataFromSamples(phenotypes, startCoord, endCoord, referenceManager, currentVariant); //@E Load in bam files
        readManager.CreateInsertionArrays(phenotypes,startCoord); // Create the insertion arrays
        
        //Adjust reference sequence for insertions. Add the insertions to this
        referenceManager.AdjustReferences(j, readManager.InsertionArrays); //@E this method adds the ins. to the reference J added
        
        //Create the Panel with Fragments using the extracted Data
        for(String type:phenotypesOrder){ // phenotype loop CIN3 etc
            //if(type.equals("Neg_Control")){
                 
                int maxReadCountForPhenotype = -1;
                
                int panelSize = referenceManager.AdjustedReferenceData.get(type).size(); //@E the number of insertions at a position
                panelManager.GetPanelFromPhenotype(type).Fragments = new HashMap[panelSize]; 
                
                Map<String,FragmentNode>[] tempFrags = new HashMap[panelSize]; //@E becomes the fragment, before the connections are made
                //@E the total number of columns is the panel size, so there are this number of fragments
                
                //int totalInsertColumns = 0;
                int panelColumnNo = 0;
                for(int columnNo = 0;columnNo<readManager.InsertionArrays.get(type).length;columnNo++){ //hmmm

                    int insertCount = readManager.InsertionArrays.get(type)[columnNo];
                    int readCountForColumn = 0;
                    //Loop through all the samples for the phenotype
                    for(int sampleNo = 0;sampleNo<phenotypes.get(type).size();sampleNo++){

                        //Extract Reads for sample
                        int sampleReadCount = readManager.GetReadsForSample(phenotypes.get(type).get(sampleNo).FileName).size();
                        readCountForColumn+=sampleReadCount;
                        //Loop through all the reads for the sample
                        for(int readNo = 0;readNo<sampleReadCount;readNo++){ // each base from bam file?? Think so, as many numbers under each point on the graph
                            gv15.Read currentRead = readManager.GetReadsForSample(
                                phenotypes.get(type).get(sampleNo).FileName).get(readNo);
                            
                            String[] readBases = currentRead.BaseValues; //@E all of the bases from each read
                            
                            //Ensure that the read is within the target region
                            int baseIndex = (startCoord - (currentRead.StartPosition+1)) + columnNo; 

                            if(baseIndex >= 0 && baseIndex < currentRead.Length){ //@E BAM sequnce is over a range. Check that the range covers the base we are currently looking for
                                String baseVal = readBases[baseIndex];
                                FragmentNode tempFragNode = new FragmentNode();  
                                    
                                if(tempFrags[panelColumnNo] == null) //@E if no fragments for this column then make a new one
                                    tempFrags[panelColumnNo] = new HashMap();
                                    
                                if(!tempFrags[panelColumnNo].containsKey(baseVal)) //@E if this base values has not already been found then put it in
                                    tempFrags[panelColumnNo].put(baseVal, tempFragNode);
                                    
                                //Increment the Read Count
                                tempFrags[panelColumnNo].get(baseVal).ReadCount++; //@E as we have found another read

                                //Get the Insert features of the current Read
                                ArrayList<InsertFeature> insertFeatures = readManager.GetInsertsForReadAtPosition(currentRead,phenotypes.get(type).get(sampleNo).FileName,columnNo,startCoord);
                                
                                //The last fragment does not have any connected fragments
                                if(panelColumnNo < tempFrags.length-1){ //@E if before the end of the column
                                if(tempFrags[panelColumnNo].get(baseVal).ConnectedFragments == null)
                                    tempFrags[panelColumnNo].get(baseVal).ConnectedFragments = new HashMap();                               
                                
                                //IF No inserts for this read therefore connect directly with the next Base
                                if(insertFeatures.isEmpty()){
                                    //Add connected Fragments
                                    if( (baseIndex+1) < currentRead.Length){
                                        String nextBaseVal = readBases[baseIndex+1];

                                        if(!tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.containsKey(nextBaseVal))
                                            tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.put(nextBaseVal, new HashSet());

                                        int connectionColumn = panelColumnNo+readManager.InsertionArrays.get(type)[columnNo]+1;
                                        if(connectionColumn < tempFrags.length)
                                            tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.get(nextBaseVal).
                                                add(connectionColumn);
                                    }                                        
                                }else{ //@E if an insertions is connected
                                    for(InsertFeature insFeature:insertFeatures){ // for as many bases are connected as insertions
                                        //Add the inserted Fragments
                                        String finalConnectedBase = null;
                                        if( (baseIndex+1) < currentRead.Length)
                                            finalConnectedBase = readBases[baseIndex+1];

                                        if(insFeature.InsertedBases.size() == 6) //@E only 5 possibles bases (why not larger than 5???)
                                            System.err.println("");
                                        
                                        AddInsertedBases(tempFrags, panelColumnNo+1, insFeature.InsertedBases, //@E call method at the bottom
                                                finalConnectedBase,panelColumnNo+readManager.InsertionArrays.get(type)[columnNo]+1);
                                        
                                        //Connect the current Fragment to the inserted Fragment
                                        String nextBaseVal = insFeature.InsertedBases.get(0);
                                       
                                        if(!tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.containsKey(nextBaseVal))
                                            tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.put(nextBaseVal, new HashSet());

                                        if( (panelColumnNo+1) < tempFrags.length) //@E only add the fragment if it does not go past the end of the graph
                                            tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.get(nextBaseVal).
                                                    add(panelColumnNo+1);
                                    }
                                }
                            }
                            }
                        }//End Read loop
                    }//End sample Loop

                    //if(insertCount > 0)
                        //System.err.println("");
                    panelColumnNo+=(insertCount+1);
                    //totalInsertColumns+=insertCount;
                    
                    if(readCountForColumn > maxReadCountForPhenotype)
                        maxReadCountForPhenotype = readCountForColumn;
                    
                }//End Column Loop    
                
                // calculate NormalizeReadCount for each FragmentNode
                double TotalReadCount, ReadCount;
                int NoEdges =0;
                for(int columnNo = 0;columnNo<panelColumnNo;columnNo++){
                    TotalReadCount = 0.0;
                    //if(columnNo == 0){
                        for(int baseType = 0;baseType<5;baseType++){
                            
                            String BaseVal = UtilityFunctions.getInstance().RowNumberToBaseType(baseType);
                            if(tempFrags[columnNo].containsKey(BaseVal)){
                                
                                TotalReadCount = TotalReadCount + tempFrags[columnNo].get(BaseVal).ReadCount;
                                
                                // find out number of edges in the 2D graph
                                for(int NextbaseType = 0;NextbaseType<5;NextbaseType++){

                                    String NextBaseVal = UtilityFunctions.getInstance().RowNumberToBaseType(NextbaseType);
                                    if(tempFrags[columnNo].containsKey(NextBaseVal)){
                                        if(tempFrags[columnNo].get(NextBaseVal).ConnectedFragments.containsKey(BaseVal)){
                                            NoEdges++;
                                        }
                                    }
                                }
                            }
                        }
                        for(int baseType = 0;baseType<5;baseType++){
                            
                            String BaseVal = UtilityFunctions.getInstance().RowNumberToBaseType(baseType);
                            if(tempFrags[columnNo].containsKey(BaseVal)){
                                
                                ReadCount = (double) tempFrags[columnNo].get(BaseVal).ReadCount;
                                tempFrags[columnNo].get(BaseVal).NormalizeReadCount = ReadCount / TotalReadCount;
                                // tempFrags[columnNo].get(BaseVal).NormalizeReadCount = (double) ReadCount;
                            }
                        }
//                    }else{
//
//                        for(int baseType = 0;baseType<5;baseType++){
//                            
//                            TotalReadCount = 0.0;                            
//                            String BaseVal = UtilityFunctions.getInstance().RowNumberToBaseType(baseType);
//                            if(tempFrags[columnNo].containsKey(BaseVal)){
//                                for(int LastcolumnNo = columnNo - 1;LastcolumnNo>=0;LastcolumnNo--){
//                                    for(int LastbaseType = 0;LastbaseType<5;LastbaseType++){
//
//                                        String LastBaseVal = UtilityFunctions.getInstance().RowNumberToBaseType(LastbaseType);
//                                        if(tempFrags[LastcolumnNo].containsKey(LastBaseVal)){
//                                            if(tempFrags[LastcolumnNo].get(LastBaseVal).ConnectedFragments.containsKey(BaseVal)){
//
//                                                if (tempFrags[LastcolumnNo].get(LastBaseVal).ConnectedFragments.get(BaseVal).contains(columnNo)) {
//                                                    TotalReadCount = TotalReadCount + tempFrags[LastcolumnNo].get(LastBaseVal).ReadCount;
//                                                }
//                                            }
//                                        }
//                                    }
//                                }
//                                ReadCount = (double) tempFrags[columnNo].get(BaseVal).ReadCount;
//                                tempFrags[columnNo].get(BaseVal).NormalizeReadCount = ReadCount / TotalReadCount;
//
//                            }
//                        }
//                    }
                }
                
                System.out.println(type + " NoEdges = " + NoEdges);
                
                if(maxReadCountForPhenotype > maxReadCount)
                    maxReadCount = maxReadCountForPhenotype;
                
                //Add the fragments to the panel fragments //@E type is the control, CIN3 etc. 
                panelManager.GetPanelFromPhenotype(type).Fragments = tempFrags;//@E fragments passed to the pannel
                panelManager.MaxReadCount = maxReadCount; //@E runs 4 times, something for each  group CIN3 etc


                int adjustedPos = 0; //@E start from the beginning 
                int addedVal = 0;
                for(int i = 0;i<flank+1;i++){ // +1 for the target base conversion identification
                    while(referenceManager.AdjustedReferenceData.get(type).get(adjustedPos).equals("INS")){
                        addedVal++; //@E number increments while the adjusted data equals INS
                        adjustedPos++;   //@E experiment with this value to see if the postion of the red line changes
                    }
                    adjustedPos++; //@E advance the position place for every base pased not just the ins
                    //@E the while is never reached again if the the adjusted Pos is not incrmented
                    
                    loopCount = loopCount + 1;
                }
                referenceManager.ShiftVals.put(type, addedVal); //@E stores the RED line value. It needs to be shifted after the insertions are added
                
            //}//End Type Check
        }//End Phenotype Loop
        //red.createFile("RED"); // print the file
    }
    
    //@E called from process fragments
    public void AddInsertedBases(Map<String,FragmentNode>[] fragments, int index,
            ArrayList<String> insertedBases,String finalConnectedBase,int finalConnectedColumn){
        
        for(int insertIndex = 0;insertIndex<insertedBases.size();insertIndex++){
        
            //Add the inserted Base
            if(fragments[index+insertIndex] == null) //@E if no fragment for this insert allready
                fragments[index+insertIndex] = new HashMap();

            FragmentNode tempFragNode = new FragmentNode(); //@E new node for the insert as it is a new column
            String insertedBase = insertedBases.get(insertIndex);
            
            if(!fragments[index+insertIndex].containsKey(insertedBase))
                fragments[index+insertIndex].put(insertedBase, tempFragNode);
                                    
            //Increment the Read Count
            fragments[index+insertIndex].get(insertedBase).ReadCount++;
            
            if(fragments[index+insertIndex].get(insertedBase).ConnectedFragments == null)
                fragments[index+insertIndex].get(insertedBase).ConnectedFragments = new HashMap();

            //Add the next Base as the connected Fragment
            if(insertIndex < insertedBases.size()-1){
                String nextInsertedBase = insertedBases.get(insertIndex+1);
                                                                       
                if(!fragments[index+insertIndex].get(insertedBase).ConnectedFragments.containsKey(nextInsertedBase))
                    fragments[index+insertIndex].get(insertedBase).ConnectedFragments.put(nextInsertedBase, new HashSet());
                
                if( (index+insertIndex+1) < fragments.length )
                    fragments[index+insertIndex].get(insertedBase).ConnectedFragments.get(nextInsertedBase).
                        add(index+insertIndex+1);
                
            }else if (insertIndex == insertedBases.size()-1){
                //Add the next non-insert base as the connected fragment of the last Insert
                
                if(!fragments[index+insertIndex].get(insertedBase).ConnectedFragments.containsKey(finalConnectedBase))
                    fragments[index+insertIndex].get(insertedBase).ConnectedFragments.put(finalConnectedBase, new HashSet());

                if( finalConnectedColumn < fragments.length ) //@E only add this if it does not go past the end
                    fragments[index+insertIndex].get(insertedBase).ConnectedFragments.get(finalConnectedBase).
                        add(finalConnectedColumn);                
            }
        }                  
    }

    public void FragmentPrinter(Map<String,FragmentNode>[] fragments){
        System.out.println("Printing Fragments\n");
        
        for(int index = 0;index<fragments.length;index++){
            if(fragments[index]!=null){
                for(String baseVal:fragments[index].keySet()){
                    System.out.print(baseVal + "(" + fragments[index].get(baseVal).ReadCount + ") ");                    
                    for(String connectedFrag:fragments[index].
                            get(baseVal).ConnectedFragments.keySet()){
                        System.out.print(connectedFrag);
                        for(Object connectedIndex:fragments[index].get(baseVal).ConnectedFragments.get(connectedFrag)){
                            System.out.print("["+connectedIndex+"]");
                        }
                    }                    
                    System.out.println("");
                }
            }
            System.out.println("");
        }
    }
}
