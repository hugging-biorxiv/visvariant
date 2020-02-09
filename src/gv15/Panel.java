package gv15;

import com.sun.javafx.geom.BaseBounds;
import com.sun.javafx.geom.transform.BaseTransform;
import com.sun.javafx.jmx.MXNodeAlgorithm;
import com.sun.javafx.jmx.MXNodeAlgorithmContext;
import com.sun.javafx.sg.prism.NGNode;
import gv15.Filters.IFilter;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import javafx.geometry.Pos;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.CubicCurve;
import javafx.scene.shape.CubicCurveTo;
import javafx.scene.shape.Line;
import javafx.scene.shape.MoveTo;
import javafx.scene.shape.Path;
import javafx.scene.shape.Shape;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
/**
 * Author Ranasi (2016)
 * Author Eisa Anwar (2017)
 */
public class Panel {
    public String PanelName;
    public Map<String,FragmentNode>[] Fragments;
    public double PositionX;
    public double PositionY;
    public int FragmentXOffset;
    public int Flank;
    private int lineThicknessMultiplier;
    private int fontMultiplier;
    private int phenotypeSize;
    
    
    private int columns;
    private int rows;
    private double columnWidth;
    private double rowHeight;
    private int renderColumns;
    private ArrayList<IFilter> panelFilters; 
            
    public Panel(String ID,double positionX,double positionY, int columns,int renderColumns,int rows,
            double columnWidth,double rowHeight, int phenotypeCount, int flank,int xOffset, int lineThicknessMulti, int fontMulti){
        this.PanelName = ID;
        this.PositionX = positionX;
        this.PositionY = positionY;
        this.columns = columns;
        this.rows = rows;
        this.columnWidth = columnWidth;
        this.rowHeight = rowHeight;
        this.phenotypeSize = phenotypeCount;
        this.Flank = flank;
        this.FragmentXOffset = xOffset;
        this.panelFilters = new ArrayList();
        this.renderColumns = renderColumns;
        this.lineThicknessMultiplier = lineThicknessMulti;
        this.fontMultiplier = fontMulti;
    }
    
    public Map<String,FragmentNode>[] getFragments(){
        return Fragments;
    }
    
    public void AddFilter(IFilter newFilter){
        panelFilters.add(newFilter);
    }
    
    public void FilterData(ArrayList<String> refereneceData) {
        //Exceute all the filters
        for(IFilter filter:panelFilters)
            filter.FilterPanel(refereneceData, Fragments);      
    }
    
    public void RenderPanel(int flank, Group renderGroup, ArrayList<String> rawReference,int maxReadCount,
            int offset, ArrayList<String> referenceData, ArrayList<Integer> requiredMoves, 
            ArrayList<Integer> insertionArray, VariantManager variantManager){
        
        //@E the rawData is not filtered
        //@E the referenceData is filtered
        
        int counter = 0;
        int noOfBlanks = 0, noOfInserts = 0;
            
        for(int i=1; i<referenceData.size(); i++) {
            if (i >= referenceData.size()) {
                break;
            }        
            if (counter >= requiredMoves.size()) {
                break;
            }
            if (referenceData.get(i).equals("INS")) {
                continue;
            }
            noOfBlanks = requiredMoves.get(counter); // number of blanks required
            for (int j=0; j<noOfBlanks; j++) {
                referenceData.add(i, "SPACE"); // add the blanks(s) after the insertions following a base
                i = i + 1;            
            }
            counter = counter + 1;
        }
         
        int totalOffset = Flank+offset;
        for(int colNum=0;colNum<totalOffset;colNum++){
            if(Fragments[colNum].isEmpty())
                offset--;
        }
        
        for(int i=0; i<flank; i++) { // offset incremented by the total number of moves to the middle
            if(i >= requiredMoves.size()) {
                break;
            }
            offset = offset + requiredMoves.get(i);
        }
        
        // find out the length of the variant
        int LengthREF = 0, LengthALT = 0;
        if(variantManager.getSelectedVariant() instanceof VariantContext){
            VariantContext tempContext = (VariantContext)variantManager.getSelectedVariant();
            for(Allele curAllele:tempContext.getAlleles()){
                // System.out.print("debug");
                if(curAllele.isReference()){
                    LengthREF = curAllele.length();
                }else{
                    if(curAllele.length() > LengthALT){
                        LengthALT = curAllele.length();
                    }
                }
            }
        }
        
        ArrayList<Shape> renderVariance;
        if(LengthREF > flank - 1 || LengthALT > flank - 1){
            renderVariance = SetupVariance(PositionX,PositionY, Flank+offset, Flank+offset+1, columnWidth, rowHeight);
        }else if(LengthREF == LengthALT && LengthREF == 1){
            renderVariance = SetupVariance(PositionX,PositionY, Flank+offset, Flank+offset+LengthREF, columnWidth, rowHeight);
        } else if (LengthREF < LengthALT || LengthREF > LengthALT || (LengthREF == LengthALT && LengthREF > 1)) { // INS or DEL
            
            noOfBlanks = 0;
            for(int i=Flank+offset; i<referenceData.size(); i++) {
                if (!referenceData.get(i).equals("INS") && !referenceData.get(i).equals("SPACE")){
                    noOfBlanks++;
                    if (noOfBlanks == LengthREF + 1) {
                        counter = i;
                        break;
                    }
                }
            }
            renderVariance = SetupVariance(PositionX,PositionY, Flank+offset, counter, columnWidth, rowHeight);
        } else {
            renderVariance = SetupVariance(PositionX,PositionY, Flank+offset, Flank+offset+1, columnWidth, rowHeight);
        }
        
        ArrayList<Shape> renderArea = SetupRenderArea((rows*2), renderColumns, columnWidth, 
                rowHeight, PositionX, PositionY);        
        ArrayList<Shape> referenceRender = SetupReferenceRender(referenceData,
                PositionX,PositionY,rowHeight,columnWidth,renderColumns);  
        ArrayList<Node> fragmentRenders = SetupFragments(rawReference,
                PositionX, PositionY, rowHeight, columnWidth, maxReadCount, referenceData);
        ArrayList<Shape> panelTitle = SetupPanelInfo(PanelName, 10, PositionY+40);
        
        renderGroup.getChildren().addAll(renderArea);
        renderGroup.getChildren().addAll(referenceRender);
        renderGroup.getChildren().addAll(renderVariance);
        renderGroup.getChildren().addAll(fragmentRenders);
        renderGroup.getChildren().addAll(panelTitle);
    }    
            
    private ArrayList<Shape> SetupVariance(double gridX, double gridY, int varianceColStart, int varianceColEnd,
            double colWidth,double rowHeight){ // @E creates the RED lines
        

        ArrayList<Shape> renderItems = new ArrayList();
                
        Line tempLine = new Line();        
        tempLine.setStartX(gridX + (varianceColStart * colWidth));
        tempLine.setStartY(gridY - rowHeight);
        tempLine.setEndX(gridX + (varianceColStart * colWidth));
        tempLine.setEndY(gridY + (5 * (rowHeight*2)));
        tempLine.setStrokeWidth(1 * lineThicknessMultiplier);
        tempLine.setStroke(Color.RED);
            
        renderItems.add(tempLine);

        Line tempLine2 = new Line();
        tempLine2.setStartX(gridX + (varianceColEnd * colWidth));
        tempLine2.setStartY(gridY - rowHeight);
        tempLine2.setEndX(gridX + (varianceColEnd * colWidth));
        tempLine2.setEndY(gridY + (5 * (rowHeight*2)));
        tempLine2.setStrokeWidth(1 * lineThicknessMultiplier);
        tempLine2.setStroke(Color.RED);
        
        renderItems.add(tempLine2);

        return renderItems;
    }
    
    
    private ArrayList<Shape> SetupRenderArea(int rows,int cols, double colWidth,double rowHeight
                                            ,double startX, double startY){
        ArrayList<Shape> renderComponents = new ArrayList<>();
        
        //Piano lines
        boolean alter = false;
        for(int row = 0;row<rows;row++){
            Line tempLine = new Line();
            tempLine.setStartX(startX + 110 + (rowHeight/2));
            tempLine.setStartY(startY + 108 + (row * rowHeight) + (rowHeight/2));
            tempLine.setEndX(startX - 110 + (cols) * colWidth - (rowHeight/2));
            tempLine.setEndY(startY + 108 + (row * rowHeight) + (rowHeight/2));
            tempLine.setStrokeWidth(rowHeight * lineThicknessMultiplier);
            
            if(alter){
                tempLine.setStroke(Color.WHITE);
                alter = false;
                renderComponents.add(tempLine);
                
                Line footerLine = new Line();
                footerLine.setStartX(startX - 100);
                footerLine.setStartY(startY + (row * rowHeight) + (rowHeight));
                footerLine.setEndX(startX );
                footerLine.setEndY(startY + (row * rowHeight) + (rowHeight));
                footerLine.setStrokeWidth(1 * lineThicknessMultiplier);    
                footerLine.setStroke(Color.GREY);               
                renderComponents.add(footerLine);
            }else{
                tempLine.setStroke(Color.GHOSTWHITE);
                alter = true;
                renderComponents.add(tempLine);
            }            
           
        }
        //Table borders
        for(int col = 0; col<=cols;col++){
            Line tempLine = new Line();
            tempLine.setStartX(startX + (col * colWidth));
            tempLine.setStartY(startY);
            tempLine.setEndX(startX + (col * colWidth));
            tempLine.setEndY(startY + (rows * rowHeight));
            tempLine.setStrokeWidth(1 * lineThicknessMultiplier);
            tempLine.setStroke(Color.LIGHTGREY);
            
            renderComponents.add(tempLine);
        }

        return renderComponents;
    }
    
    
    private ArrayList<Shape> SetupPanelInfo(String panelName,double startX,
            double startY){
        ArrayList<Shape> renderElements = new ArrayList<Shape>();         
        
        Text tempTextTitle = new Text(startX, startY + 50, panelName);
        tempTextTitle.setFont(Font.font("Verdana", FontWeight.BOLD, 25 * fontMultiplier));
        tempTextTitle.setFill(Color.LIGHTCORAL);
        renderElements.add(tempTextTitle);    
        
        Text tempTextReadCount = new Text(startX, startY + 250, "#Sample: " + phenotypeSize);
        tempTextReadCount.setFont(Font.font("Verdana", FontWeight.BOLD, 20 * fontMultiplier));
        tempTextReadCount.setFill(Color.LIGHTCORAL);
        renderElements.add(tempTextReadCount);    
        
        return renderElements;
    }
    
    
    private ArrayList<Shape> SetupReferenceRender(ArrayList<String> referenceData,
            double startX,double startY,double rowHeight,double colWidth,int cols){
        ArrayList<Shape> renderTexts = new ArrayList<Shape>();
        String colText;
                
        //Reference BasePairs
        for(int refIndex = 0;refIndex<cols;refIndex++){
            
            if(refIndex < referenceData.size()){
                if (referenceData.get(refIndex).equals("SPACE")) {
                    colText = "SPC"; // Label for the blank column
                } else {               
                    colText = referenceData.get(refIndex);
                }
                // print the base in the middle
                Text tempText = new Text(startX - 10 + (colWidth*refIndex) + (colWidth/2), 
                        startY, colText); //@E base of the column

                tempText.setFont(Font.font ("Verdana", 10 * fontMultiplier));
                renderTexts.add(tempText);
            }
        }        
        
        //Read BasePair Types
        for(int baseType = 0;baseType<5;baseType++){
            Text tempText = new Text(startX - 95,
                    startY + (baseType*rowHeight*2) + 110,UtilityFunctions.
                            getInstance().RowNumberToBaseType(baseType));
            tempText.setFont(Font.font("Verdana", FontWeight.BOLD, 25 * fontMultiplier));
            renderTexts.add(tempText); //@E letters of bases on the left
        }
        
        return renderTexts;
    }
    
    private ArrayList<Node> SetupFragments(ArrayList<String> referenceData,
            double gridX, double gridY, double rowHeight, 
            double colWidth,int maxReadCount, ArrayList<String> columnReferenceData){
        
        ArrayList<Node> renderElements = new ArrayList<>();
        
        int XOFFSET = FragmentXOffset;
        int YOFFSET = 16;
        int skippedFragments = 0;   //@E fragments below the read count threshold
        int blankFragments = 0; //@E total of all blank spaces added to the sequences
        int lastPrintcol = 0;

        for(int colNum=0; colNum < renderColumns+skippedFragments; colNum++){
            if(colNum >= Fragments.length)
                break;
            
            if(Fragments[colNum].isEmpty()){
                skippedFragments++; //@E if empty dont draw on that column
                continue;
            }
            
            if(columnReferenceData.get(colNum - skippedFragments + blankFragments).equals("SPACE")) { // if continues here and does not add the blank spaces
                blankFragments = blankFragments + 1;
                colNum = colNum - 1; // there are no blanks in the fragment so they should not be skipped
                continue; //@E move fowards but don't skip a fragment, leave a blank space               
            }
                        
            for(int baseType = 0;baseType<5;baseType++){ // loop through the base possibilities
                
                if(Fragments[colNum]!=null && Fragments[colNum].containsKey(UtilityFunctions.getInstance().RowNumberToBaseType(baseType))
                        && (colNum-skippedFragments) < renderColumns){
                    
                    FragmentNode val = Fragments[colNum].get(UtilityFunctions.getInstance().RowNumberToBaseType(baseType));
                    
                    float readSize = 1 + (val.ReadCount/(maxReadCount*1.0f)*13.0f);                   
                    
                    //@E 1 fragment lines
                    Line tempLine = new Line();
                    tempLine.setStartX(( (colNum-skippedFragments) * colWidth) + (blankFragments * colWidth) + gridX + XOFFSET + 4*(readSize/2));
                    tempLine.setStartY((baseType * rowHeight * 2) + gridY + YOFFSET + 4*(readSize/2));
                    tempLine.setEndX(( (colNum-skippedFragments) * colWidth) + (blankFragments * colWidth) + colWidth + gridX - XOFFSET - 4*(readSize/2)); //@E  + skipColumn
                    tempLine.setEndY((baseType * rowHeight * 2) + gridY + YOFFSET + 4*(readSize/2));
                    tempLine.setStrokeWidth(readSize * lineThicknessMultiplier);
                   
                    if(!referenceData.get(colNum).equals("INS")){
                        if(referenceData.get(colNum).equals(UtilityFunctions.
                            getInstance().RowNumberToBaseType(baseType)))
                            tempLine.setStroke(Color.web(UtilityFunctions.getInstance().ReadColour_Unvaried));
                        else
                            tempLine.setStroke(Color.web(UtilityFunctions.getInstance().ReadColour_Varied));
                    }else
                        tempLine.setStroke(Color.web(UtilityFunctions.getInstance().ReadColour_Insertion));
                    
                    renderElements.add(tempLine);                         

                    //Connect the fragments
                    if(colNum < Fragments.length){
                        for(int nextBase = 0;nextBase<5;nextBase++){

                            if(val.ConnectedFragments != null && val.ConnectedFragments.containsKey(UtilityFunctions.getInstance().RowNumberToBaseType(nextBase))){
                                
                                String connectedVal = UtilityFunctions.getInstance().RowNumberToBaseType(nextBase);
                                HashSet<Integer> connectedColumns = val.ConnectedFragments.get(UtilityFunctions.getInstance().RowNumberToBaseType(nextBase));

                                for (Integer colVal : connectedColumns) {
                                    
                                    if(colVal >= Fragments.length)
                                        System.err.println("panel.java line = 337");
                                    
                                    if(!Fragments[colVal].containsKey(connectedVal))
                                        System.err.println("");
                                    
                                    if(Fragments[colVal].containsKey(connectedVal) &&
                                            (colVal-skippedFragments) < this.renderColumns){

                                        int connectionEndColumn = GetConnectionEnd(colNum,colVal)-skippedFragments;
                                    
                                        // if the connected column falls on a blank space then increment it
                                        while(columnReferenceData.get(connectionEndColumn + blankFragments).equals("SPACE")) {
                                            connectionEndColumn = connectionEndColumn + 1;
                                        }
                                    
                                        float nextReadSize = 1 + ((Fragments[colVal].get(connectedVal).ReadCount
                                            /(maxReadCount*1.0f))*13.0f);    
                                                                    
                                        if(nextBase == baseType){ //@E if the base has not changed then straight line
                                            Line connectorLine = new Line();
                                            //skip additional fragments
                                            connectorLine.setStartX(( (colNum-skippedFragments) * colWidth) + (blankFragments * colWidth) 
                                                    +colWidth + gridX - XOFFSET - 4*(readSize/2));
                                            connectorLine.setEndX(( (connectionEndColumn) * colWidth) + (blankFragments * colWidth)
                                                    + gridX + XOFFSET + 4*(readSize/2));
                                            
                                            if(nextReadSize<readSize){                                      
                                                connectorLine.setStrokeWidth(nextReadSize * lineThicknessMultiplier);
                                                connectorLine.setStartY((baseType * rowHeight * 2) + gridY + YOFFSET + 4*(nextReadSize/2));
                                                connectorLine.setEndY((baseType * rowHeight * 2) + gridY + YOFFSET + 4*(nextReadSize/2));
                                            }else{
                                                connectorLine.setStrokeWidth(readSize * lineThicknessMultiplier);
                                                connectorLine.setStartY((baseType * rowHeight * 2) + gridY + YOFFSET + 4*(readSize/2)); //start and end
                                                connectorLine.setEndY((baseType * rowHeight * 2) + gridY + YOFFSET + 4*(readSize/2));
                                            }
                                            if(referenceData.get(colNum).equals("INS") || //@E the insertion could make the base the same
                                                    referenceData.get(colVal).equals("INS"))                                        
                                                connectorLine.setStroke(Color.web(UtilityFunctions.getInstance().ReadColour_Insertion)); // @E purple line
                                            else{
                                                if(referenceData.get(colVal).equals(UtilityFunctions.
                                                        getInstance().RowNumberToBaseType(baseType)))
                                                    connectorLine.setStroke(Color.web(UtilityFunctions.getInstance().ReadColour_Unvaried)); //@E grey line
                                                else
                                                    connectorLine.setStroke(Color.web(UtilityFunctions.getInstance().ReadColour_Varied)); // @E orange line
                                            }
                                            
                                            renderElements.add(connectorLine);  
                                        }else{ //@E there needs to be a curved line as the base has changed
                                            //Join the fragments
                                            CubicCurve tempCurve = new CubicCurve(
                                                ( (colNum-skippedFragments) * colWidth) + (blankFragments * colWidth) + colWidth + gridX - XOFFSET,//@E start X point
                                                (baseType * rowHeight * 2) + gridY + YOFFSET +1,// - 2*(readSize/2), //@E start Y coordinate
                                                    
                                                ( (colNum-skippedFragments) * colWidth) + (blankFragments * colWidth) + colWidth + gridX - XOFFSET + 72,
                                                (baseType * rowHeight * 2) + gridY + YOFFSET +1,// - 2*(readSize/2), //@E +1 to Y to line up lines

                                                ( (connectionEndColumn) * colWidth) + (blankFragments * colWidth) + gridX + XOFFSET  - 72,
                                                (nextBase * rowHeight * 2) + gridY + YOFFSET+1,// + 2*(nextReadSize/2),

                                                ( (connectionEndColumn) * colWidth) + (blankFragments * colWidth) + gridX + XOFFSET, //@E end X  + skipColumn
                                                (nextBase * rowHeight * 2) + gridY + YOFFSET+1// + 2*(nextReadSize/2)
                                            );
                                            if(referenceData.get(colNum).equals("INS") 
                                                    || referenceData.get(colVal).equals("INS"))
                                                tempCurve.setStroke(Color.web(UtilityFunctions.getInstance().ReadColour_Insertion));
                                            else{
                                                if(referenceData.get(colNum).equals(UtilityFunctions.
                                                        getInstance().RowNumberToBaseType(baseType)) &&
                                                        referenceData.get(colVal).equals(UtilityFunctions.
                                                        getInstance().RowNumberToBaseType(nextBase)))
                                                    tempCurve.setStroke(Color.web(UtilityFunctions.getInstance().ReadColour_Unvaried));
                                                else
                                                    tempCurve.setStroke(Color.web(UtilityFunctions.getInstance().ReadColour_Varied));
                                            }
                                            tempCurve.setFill(null);
                                            renderElements.add(tempCurve);      
                                            //@E makes the curve thicker at the fragment nodes so the cubic curves thickness is the same as the thickness of the fragment lines
                                            //readsize times by the linethickness multiplier
                                            Path path = new Path();
                                            for (int i = 0; i <= (nextReadSize*lineThicknessMultiplier)-1; i++) {
                                                path.getElements().addAll(
                                                        new MoveTo(tempCurve.getStartX(), tempCurve.getStartY()),
                                                        new CubicCurveTo(
                                                                tempCurve.getControlX1() ,
                                                                tempCurve.getControlY1() ,
                                                                tempCurve.getControlX2() ,
                                                                tempCurve.getControlY2() ,
                                                                tempCurve.getEndX() ,
                                                                tempCurve.getEndY() + i 
                                                        )
                                                );
                                            }
                                            path.setStroke(tempCurve.getStroke());
                                            
                                            Path path2 = new Path();
                                            for (int i = 0; i <= (readSize*lineThicknessMultiplier)-1; i++) {
                                                path2.getElements().addAll(
                                                        new MoveTo(tempCurve.getStartX(), tempCurve.getStartY() + i),
                                                        new CubicCurveTo(
                                                                tempCurve.getControlX1() ,
                                                                tempCurve.getControlY1() ,
                                                                tempCurve.getControlX2() ,
                                                                tempCurve.getControlY2() ,
                                                                tempCurve.getEndX() ,
                                                                tempCurve.getEndY() 
                                                        )
                                                );
                                            }
                                            path2.setStroke(tempCurve.getStroke());

                                            renderElements.add(path);
                                            renderElements.add(path2);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    //Add Read count render
                    // Text tempText = new Text(0,0,Integer.toString(val.ReadCount )); //@E number of fragments
                    Text tempText = new Text(0,0,String.format("%.3f", val.NormalizeReadCount)); //@E number of fragments
                    double width = 4*(tempText.getLayoutBounds().getWidth());
                    tempText.setFont(Font.font ("Verdana", 9 * fontMultiplier));
                    VBox textBox = new VBox();
                    textBox.getChildren().addAll(tempText);
                    textBox.setAlignment(Pos.BASELINE_CENTER);
                    textBox.setLayoutX(( (colNum-skippedFragments) * colWidth) + (blankFragments * colWidth) + gridX + 10 + (colWidth/2) - (width/2));
                    textBox.setLayoutY((baseType * rowHeight * 2) + gridY + YOFFSET + 7*(readSize/2));                    
                    renderElements.add(textBox);  
                    
                }
            }
        }
        
        return renderElements;
    }
    
    private int GetConnectionEnd(int currentColumn,int targetColumn){  //@E the number of columns the connection is spread over
        int endColumn = targetColumn;
        
        while(currentColumn<targetColumn-1){
            if(Fragments[currentColumn+1].isEmpty()) //@E lok through the fragments to find where the next one is
                endColumn--;
            
            currentColumn++;
        }

        return endColumn;
    }
    
}



