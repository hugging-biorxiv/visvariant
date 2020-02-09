/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gv15;

import gv15.Filters.ReadCountFilter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import javafx.scene.Group;

/**
 * Author Ranasi (2016)
 * Author Eisa Anwar (2017)
 */
public class PanelManager {
    
    private ArrayList<Panel> enginePanels;
    private Map<Panel,ArrayList<String>> referenceData = new HashMap<Panel,ArrayList<String>>(); //@E Map of the filtered data
    private Map<Panel,ArrayList<Integer>> insertionArrays = new HashMap<Panel,ArrayList<Integer>>(); // Map to store the number of insertions with each base for each panel
    private Map<Panel,ArrayList<Integer>> requiredMovesMap = new HashMap<Panel,ArrayList<Integer>>(); // Map to tore the required moves by each colum of each panel
    public int MaxReadCount;
      
    public PanelManager(){
        enginePanels = new ArrayList();
        MaxReadCount = 0;
    }
    
    public void AddPanel(String panelName,double startX,double startY, int phenotypeSize, int flank,
            double columnWidth,double rowHeight,int xOffset,int renderColumns, int lineMultiplier, int fontMultiplier){
            
        Panel tempPanel = new Panel(panelName, startX, startY, (flank*2)+1,renderColumns,5,
            columnWidth,rowHeight, phenotypeSize, flank,xOffset, lineMultiplier, fontMultiplier);
        
        //Add filters
        tempPanel.AddFilter(new ReadCountFilter());
        
        enginePanels.add(tempPanel);
    }    
    
    public void RenderPanels(int flank, Group root,ReferenceManager referenceManager, VariantManager variantManager){
        for(Panel panel:enginePanels){
            //if(panel.PanelName.equals("Neg_Control"))
            panel.RenderPanel(flank, root, referenceManager.GetReferenceForType(panel.PanelName),
                    MaxReadCount,referenceManager.ShiftVals.get(panel.PanelName), referenceData.get(panel), 
                    requiredMovesMap.get(panel), insertionArrays.get(panel), variantManager); //@E raw referencs added and sent
        }    
    }    
    
    public Panel GetPanelFromPhenotype(String type){
        for(Panel panel:enginePanels){
            if(panel.PanelName.equals(type))
                return panel;
        }       
        
        return null;
    }    
    
    public void RunFilters(ReferenceManager referenceManager) {
        for(Panel panel:enginePanels){ // for each panel           
            ArrayList<String> rawReference = new ArrayList(referenceManager.GetReferenceForType(panel.PanelName));
            ArrayList<String> reference = new ArrayList<String>();
            
            for(String value : rawReference) { // copy the values
                reference.add(value);                
            }         
            
            panel.FilterData(reference); // call method to run filters            
            referenceData.put(panel, reference);
        }       
    }
    
    public void AlignPanels(int flank, HashMap<String,ArrayList<String>> AdjustedReferenceData, HashMap<String,ArrayList<Phenotype>> phenotypes, int colRender) { //@E line up all the fragments between the panels
       
        ArrayList<String> sequence;
        ArrayList<Integer> insertionArray; 
        int insCount = 0;
        
        for(Panel panel:enginePanels) {           
            sequence = referenceData.get(panel);
            insertionArray = new ArrayList<Integer>(); // new arraylist each time to prevent overwriting
            for (int i=1; i<sequence.size(); i++) {              
                while(sequence.get(i).equals("INS")) { // find how many insertions after this
                    insCount = insCount + 1;
                    if (i < (sequence.size()-1)) { // if there are still more bases in the sequence
                        i = i + 1; // increment the main counter as this base was an insertion
                    } else { // there are no more bases in the sequence
                        break; // leave the while loop, don't check any more bases even
                    }
                }
                insertionArray.add((Integer)insCount);  // how to put values in the same spot each time , i is skipped
                insCount = 0; // reset the counter                  
            }
            insertionArrays.put(panel, insertionArray); // put the array in the map            
        }        
        
        ArrayList<Integer> largestInsertion = new ArrayList<Integer>(); // largest number of insertions for a panel for eeach base
               
        for (int i=0; i<=flank*2; i++) { // calculare the largest insertion of each base in the variant        
            largestInsertion.add(i,0); // add arbitary value so a comparison can be made
            for(Panel panel:enginePanels) {
                if (i >= insertionArrays.get(panel).size()) {
                    break;
                }                
                if (largestInsertion.get(i) < insertionArrays.get(panel).get(i)) {
                    largestInsertion.remove(i); // remove the existing value otherwise it remains in the list
                    largestInsertion.add(i, insertionArrays.get(panel).get(i)); // put the larger value in
                }
            }
        }
        
        ArrayList<Integer> requiredMoves;
        int requiredMove;
        for(Panel panel:enginePanels) {
            requiredMoves = new ArrayList<Integer>();
            
            for(int i=0; i<flank*2; i++) {
            if (i >= insertionArrays.get(panel).size()) {
                break;
            }
                requiredMove = largestInsertion.get(i) - insertionArrays.get(panel).get(i); // calculate the moves required to align
                requiredMoves.add(requiredMove); // will not be smaller than 0
            }
            requiredMovesMap.put(panel, requiredMoves);
        }
    }  
}