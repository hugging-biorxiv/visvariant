package gv15;

import com.opencsv.CSVReader;
//import data.cache.ConsensusFileCache;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.application.Application;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.SnapshotParameters;
import javafx.scene.SnapshotResult;
import javafx.scene.image.WritableImage;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.stage.Stage;
import javafx.util.Callback;
import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.stream.FileImageOutputStream;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

/**
 * Author Ranasi (2016)
 * Author Eisa Anwar (2017 - 2018)
 */
public class Engine{
    
    //Engine Settings
    public double WIDTH, HEIGHT;
    public int FLANK;
    public String DataPath;
    public String CachePath;
    public String OutputPath;
    public String ReferencePath;
    public String VariantPath;
    public String PhenotypePath;
    public String OutputType;
    public String TestFilePath;

    public double GridStartX;
    public double GridStartY;
    public double PanelSeparation;
    public int RenderColumns;
    public int ColumnWidth;
    public int RowHeight;
    public int FragmentXOffset;  
    public int PhenotypeColumn;   
    public int fontMultiplier;
    public int lineMultiplier;
    //public int minimumVariantNumber; //@E minimun number of variants needed to render the fragment
    
    //Debug Commands
    boolean TESTINGPANELS = false;
    boolean TESTING = false;
    
    //Components
    DataManager dataManager;
    PanelManager panelManager[];
    VariantManager variantManager;
    FragmentManager fragmentManager;
    ReferenceManager referenceManager[];
    //TabletDataHandler tabletDataHandler;
    HashMap<String,ArrayList<Phenotype>> phenotypes = new HashMap();
    ArrayList<String> phenotypesOrder = new ArrayList();
    
    public Engine(String[] args){

        if(!TESTINGPANELS){
            ParseArguments(args);
            
            System.out.println("Starting Engine");
            //Setup Import Utils
            
            dataManager = new DataManager(DataPath,VariantPath,PhenotypePath,PhenotypeColumn);
            dataManager.ImportPhenotypes(phenotypes);
            dataManager.ImportPhenotypesOrder(phenotypesOrder);
            
            //Setup Variants
            variantManager = new VariantManager(dataManager.ImportVCFFile());     
            
            panelManager = new PanelManager[variantManager.TotalVariantCount];
            referenceManager = new ReferenceManager[variantManager.TotalVariantCount];
            
            //Loop through all the Variants
            for(int i = 0;i < variantManager.TotalVariantCount;i++){
                
                variantManager.setVariant(i);
                System.out.println("Processing variant "+variantManager.getSelectedVariant().getStart());
                
                //Setup Panels
                panelManager[i] = new PanelManager();
                referenceManager[i] = new ReferenceManager(ReferencePath); //@E array of this object

                int count = 0;
                for(String type:phenotypesOrder){ // the 4 current panels
                    panelManager[i].AddPanel(type, GridStartX, GridStartY + (PanelSeparation*count), phenotypes.get(type).size(),
                            FLANK, ColumnWidth, RowHeight,FragmentXOffset,RenderColumns, lineMultiplier, fontMultiplier);  //@E create the image
                    count++;
                }

                // read in the human reference genome
                
                //Setup Fragments
                fragmentManager = new FragmentManager(DataPath);               
                try {
                    fragmentManager.ProcessFragments(i, phenotypes, phenotypesOrder, referenceManager[i], //@E i added. This is done 5 times 
                            panelManager[i],FLANK,variantManager.getSelectedVariant());
                } catch (Exception ex) {
                    Logger.getLogger(Engine.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            
        }
    }

    public void Render(Stage stage){
        System.out.println("Rendering");
        //if(TESTINGPANELS)
            //CreateTestPanel();
        
        //Render for all Variants
        for(int i = 0;i<variantManager.TotalVariantCount;i++){ //@E for all the tables
            variantManager.setVariant(i);
            Group root = new Group();
            
            panelManager[i].RunFilters(referenceManager[i]);

            panelManager[i].AlignPanels(FLANK, referenceManager[i].AdjustedReferenceData, phenotypes, RenderColumns);

            panelManager[i].RenderPanels(FLANK, root, referenceManager[i], variantManager);
           
            // maximum offset
            if(!TESTINGPANELS)
                root.getChildren().add(SetupChartTitle(20, 100));

            Scene scene = new Scene(
                    root,
                    WIDTH, HEIGHT,
                    Color.rgb(255,255,255)
            );

            stage.setScene(scene);
            stage.setResizable(true);
            //stage.show(); //@E causes the fisrt image to not be fully rendered. (size only of screen resolution)

            OutputResultsToImage(scene);
        }
        
        //stage.setMaximized(true);        
    }
    
    private void ParseArguments(String[] args){
        SetupDefaultParameters();
      
        //No arguments entered therefore show OPTIONS and exit
        if(args.length == 0){
            System.out.println(UtilityFunctions.getInstance().AppOptions.OptionsList);
            if (!TESTING) System.exit(0);
        }
        if (!TESTING) {
            for(int argNum = 0;argNum<args.length;argNum++){
                int endIndex = args[argNum].indexOf("=");
                int startIndex = args[argNum].indexOf("-");

                String argType = args[argNum].substring(startIndex+1,endIndex);
                String parameterVal = args[argNum].substring(endIndex+1);

                //If a prefs path is provided extract settings from the target file
                if(argType.equals("pref") || argType.equals("-pref")){
                    SetPrefsFile(parameterVal);
                    break;
                }
                // System.out.println(parameterVal);

                //Handle options
                if(argType.equals("datapath") || argType.equals("-datapath")){
                    DataPath = parameterVal;
                    System.out.println("datapath = " + DataPath);
                }else if(argType.equals("referencepath") || argType.equals("-referencepath")){
                    ReferencePath = parameterVal;
                    System.out.println("referencepath = " + ReferencePath);
                }else if(argType.equals("variantpath") || argType.equals("-variantpath")){
                    VariantPath = parameterVal;
                    System.out.println("variantpath = " + VariantPath);
                }else if(argType.equals("phenotypepath") || argType.equals("-phenotypepath")){
                    PhenotypePath = parameterVal;
                    System.out.println("phenotypepath = " + PhenotypePath);
                }else if(argType.equals("outputpath") || argType.equals("-outputpath")){
                    OutputPath = parameterVal;
                    System.out.println("outputpath = " + OutputPath);
                }else if(argType.equals("r") || argType.equals("-r")){
                    UtilityFunctions.getInstance().VariantCoordinate = parameterVal;
                    System.out.println("r = " + UtilityFunctions.getInstance().VariantCoordinate);
                }else if(argType.equals("width") || argType.equals("-width")){
                    WIDTH = Double.parseDouble(parameterVal);
                    System.out.println("width = " + WIDTH);
                }else if(argType.equals("height") || argType.equals("-height")){
                    HEIGHT = Double.parseDouble(parameterVal);
                    System.out.println("height = " + HEIGHT);
                }else if(argType.equals("flank") || argType.equals("-flank")){
                    FLANK = Integer.parseInt(parameterVal); 
                    System.out.println("flank = " + FLANK);
                }else if(argType.equals("gridstartx") || argType.equals("-gridstartx")){
                    GridStartX = Double.parseDouble(parameterVal);
                    System.out.println("gridstartx = " + GridStartX);
                }else if(argType.equals("gridstarty") || argType.equals("-gridstarty")){
                    GridStartY = Double.parseDouble(parameterVal);
                    System.out.println("gridstarty = " + GridStartY);
                }else if(argType.equals("panelseparation") || argType.equals("-panelseparation")){
                    PanelSeparation = Double.parseDouble(parameterVal);
                    System.out.println("panelseparation = " + PanelSeparation);
                }else if(argType.equals("columns") || argType.equals("-columns")){
                    RenderColumns = Integer.parseInt(parameterVal);
                    System.out.println("columns = " + RenderColumns);
                }else if(argType.equals("columnwidth") || argType.equals("-columnwidth")){
                    ColumnWidth = Integer.parseInt(parameterVal);
                    System.out.println("columnwidth = " + ColumnWidth);
                }else if(argType.equals("rowheight") || argType.equals("-rowheight")){
                    RowHeight = Integer.parseInt(parameterVal);
                    System.out.println("rowheight = " + RowHeight);
                }else if(argType.equals("fragmentxoffset") || argType.equals("-fragmentxoffset")){
                    FragmentXOffset = Integer.parseInt(parameterVal);
                    System.out.println("fragmentxoffset = " + FragmentXOffset);
                }else if(argType.equals("outputtype") || argType.equals("-outputtype")){
                    OutputType = parameterVal;
                    System.out.println("outputtype = " + OutputType);
                }else if(argType.equals("phenotypecolumn") || argType.equals("-phenotypecolumn")){
                    PhenotypeColumn = (Integer.parseInt(parameterVal)-1);
                    System.out.println("phenotypecolumn = " + PhenotypeColumn);
                }else if(argType.equals("fontSizeMultiplier") || argType.equals("-fontSizeMultiplier")){
                    fontMultiplier = Integer.parseInt(parameterVal);
                    System.out.println("fontSizeMultiplier = " + fontMultiplier);
                }else if(argType.equals("lineThicknessMultiplier") || argType.equals("-lineThicknessMultiplier")){
                    lineMultiplier = Integer.parseInt(parameterVal);
                    System.out.println("lineThicknessMultiplier = " + lineMultiplier);
                }else if(argType.equals("readcountrenderthreshold") || argType.equals("-readcountrenderthreshold")){
                    UtilityFunctions.getInstance().ReadCountRenderThreshold = Integer.parseInt(parameterVal);
                    System.out.println("readcountrenderthreshold = " + UtilityFunctions.getInstance().ReadCountRenderThreshold);
                }else if(argType.equals("insertionsonlyatvariant") || argType.equals("-insertionsonlyatvariant")){
                    if(parameterVal.equals("0"))
                        UtilityFunctions.getInstance().InsertionsOnlyAtVariant = false;
                    else
                        UtilityFunctions.getInstance().InsertionsOnlyAtVariant = true;
                    System.out.println("insertionsonlyatvariant = " + UtilityFunctions.getInstance().InsertionsOnlyAtVariant);
                }else if(argType.equals("readcolour_unvaried") || argType.equals("-readcolour_unvaried")){
                    UtilityFunctions.getInstance().ReadColour_Unvaried = parameterVal;
                    System.out.println("readcolour_unvaried = " + UtilityFunctions.getInstance().ReadColour_Unvaried);
                }else if(argType.equals("readcolour_varied") || argType.equals("-readcolour_varied")){
                    UtilityFunctions.getInstance().ReadColour_Varied = parameterVal;
                    System.out.println("readcolour_varied = " + UtilityFunctions.getInstance().ReadColour_Varied);
                }else if(argType.equals("readcolour_insertion") || argType.equals("-readcolour_insertion")){
                    UtilityFunctions.getInstance().ReadColour_Insertion = parameterVal;
                    System.out.println("readcolour_insertion = " + UtilityFunctions.getInstance().ReadColour_Insertion);
                }
            }
        }
        // Ensure the .scri-bioinf folder exists
	// File fldr = new File(".scri-bioinf");
	// fldr.mkdirs();  

	// Cached reference file
	// ConsensusFileCache.setIndexFile(new File(fldr, "tablet-refs.xml"));
    }
    
    private void SetupDefaultParameters(){
        String BasePath;
        if (TESTING) {
            BasePath = "C:\\Users\\lau01\\Documents\\NetBeansProjects\\";
            DataPath = BasePath + "data\\test_region\\bam_filtered\\";
            OutputPath = BasePath;
            ReferencePath = BasePath + "data\\ref.fasta";
            VariantPath = BasePath + "data\\merge_bcftools_recall_fixed_subset.vcf";
            PhenotypePath = BasePath + "data\\ArtisticPhenotype.csv";
        }
        FLANK = 6;
        GridStartX = 900;
        GridStartY = 400;
        PanelSeparation = 900;
        RenderColumns = 40;
        ColumnWidth = 180;
        RowHeight = 72;
        FragmentXOffset = 60;
        WIDTH = 7680;
        HEIGHT = 4320;

        OutputType = "png";
        PhenotypeColumn = 3;
        fontMultiplier = 4;
        lineMultiplier = 4;
        UtilityFunctions.getInstance().ReadCountRenderThreshold = 50;
        UtilityFunctions.getInstance().VariantCoordinate = "chr1;871328;A;AG";
        UtilityFunctions.getInstance().InsertionsOnlyAtVariant = false;
        UtilityFunctions.getInstance().ReadColour_Unvaried = "#30dd47";
        UtilityFunctions.getInstance().ReadColour_Varied = "#23cbff";
        UtilityFunctions.getInstance().ReadColour_Insertion = "#b238ff";
    }
    
//    private static ImageView createScaledView(Node node, int scale) {
//        final Bounds bounds = node.getLayoutBounds();
//
//        final WritableImage image = new WritableImage(
//            (int) Math.round(bounds.getWidth() * scale),
//            (int) Math.round(bounds.getHeight() * scale));
//
//        final SnapshotParameters spa = new SnapshotParameters();
//        spa.setTransform(javafx.scene.transform.Transform.scale(scale, scale));
//
//        final ImageView view = new ImageView(node.snapshot(spa, image));
//        view.setFitWidth(bounds.getWidth());
//        view.setFitHeight(bounds.getHeight());
//        
//        return view;
//    }
    
    
    private void OutputResultsToImage(Scene scene){
        double scale = 0.5;
        WritableImage writableImage = new WritableImage((int) Math.round(scene.getWidth() * scale), (int) Math.round(scene.getHeight() * scale));

        scene.snapshot(writableImage);
        
        String fileName = variantManager.getSelectedVariant().getContig() + "_" + 
                Integer.toString(variantManager.getSelectedVariant().getStart()); 
        if(OutputType.equals("png")){
            BufferedImage tempImg = SwingFXUtils.fromFXImage(writableImage, null);
            File outputfile = new File(OutputPath+fileName+"_results.png");
            try{
                ImageIO.write(tempImg, "png", outputfile);
            }catch(Exception e){

            }
        }else if(OutputType.equals("jpeg")){
           File fa = new File(OutputPath+fileName+"_results.jpg");
           RenderedImage renderedImage = SwingFXUtils.fromFXImage(writableImage, null);
           BufferedImage image2 = new BufferedImage((int)WIDTH, (int)HEIGHT, BufferedImage.TYPE_INT_RGB); 
           image2.setData(renderedImage.getData());
           try{
               ImageWriter writer = (ImageWriter)ImageIO.getImageWritersByFormatName("jpeg").next();
               ImageWriteParam iwp = writer.getDefaultWriteParam();
               iwp.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
               iwp.setCompressionQuality(1);   // a float between 0 and 1
               // 1 specifies minimum compression and maximum quality
               FileImageOutputStream output = new FileImageOutputStream(fa);
               writer.setOutput(output);
               IIOImage iioimage = new IIOImage(image2, null, null);
               writer.write(null, iioimage, iwp);
               writer.dispose();

           }catch (Exception e){

           }     
        }        
    }

    private Text SetupChartTitle(int startX,int startY){
        Text details = new Text(startX,startY,"");
        details.setFont(Font.font("Verdana", FontWeight.BOLD, 25 * fontMultiplier));
        details.setFill(Color.LIGHTSLATEGREY);
        //details.setA
        
        String displayText = variantManager.getSelectedVariant().getContig() +
                ":" + variantManager.getSelectedVariant().getStart() + " ";
        
        if(variantManager.getSelectedVariant() instanceof VariantContext){
            VariantContext tempContext = (VariantContext)variantManager.getSelectedVariant();
            for(Allele curAllele:tempContext.getAlleles()){
                for(byte base:curAllele.getBases()){
                    displayText+=Character.toString ((char) base);
                }
                displayText+=">";
            }
        }

        details.setText(displayText.substring(0, displayText.length()-1));

        return details;
    }   
    
    private File SetPrefsFile(String filePath){
        //Read Engine Preferences
	try (BufferedReader br = new BufferedReader(new FileReader(filePath)))
        {
            String sCurrentLine;
            while ((sCurrentLine = br.readLine()) != null) {
                String parameterName = sCurrentLine.substring(0, sCurrentLine.indexOf("="));
                String parameterVal = sCurrentLine.substring(sCurrentLine.indexOf("=")+1);
                parameterVal = parameterVal.trim();
                
                //Setting the parameters as specified in the prefs file 
                switch(parameterName){
                    case "datapath": DataPath = parameterVal;
                        break;
                    case "cachepath": CachePath = parameterVal;
                        break;
                    case "outputpath": OutputPath = parameterVal;
                        break;         
                    case "referencepath": ReferencePath = parameterVal;
                        break;
                    case "variantpath": VariantPath = parameterVal;
                        break;     
                    case "phenotypepath": PhenotypePath = parameterVal;
                        break;
                    case "width": WIDTH = Double.parseDouble(parameterVal);
                        break;
                    case "height": HEIGHT = Double.parseDouble(parameterVal);
                        break; 
                    case "flank": FLANK = Integer.parseInt(parameterVal);
                        break; 
                    case "gridstartx": GridStartX = Double.parseDouble(parameterVal);
                        break;
                    case "gridstarty": GridStartY = Double.parseDouble(parameterVal);
                        break;
                    case "panelseparation": PanelSeparation = Double.parseDouble(parameterVal);
                        break;
                    case "columns": RenderColumns = Integer.parseInt(parameterVal);
                        break;
                    case "columnwidth": ColumnWidth = Integer.parseInt(parameterVal);
                        break;
                    case "rowheight": RowHeight = Integer.parseInt(parameterVal);
                        break;
                    case "fragmentxoffset": FragmentXOffset = Integer.parseInt(parameterVal);
                        break;
                    case "outputtype": OutputType = parameterVal;
                        break;
                    case "phenotypecolumn": PhenotypeColumn = (Integer.parseInt(parameterVal)-1);
                        break;
                    case "readcountrenderthreshold": UtilityFunctions.getInstance().ReadCountRenderThreshold 
                                                        =  Integer.parseInt(parameterVal);
                        break; 
                    case "insertionsonlyatvariant": 
                        if(parameterVal.equals("0"))
                            UtilityFunctions.getInstance().InsertionsOnlyAtVariant = false;
                        else
                            UtilityFunctions.getInstance().InsertionsOnlyAtVariant = true;
                        break;
                    case "lineThicknessMultiplier":  lineMultiplier = (Integer.parseInt(parameterVal));
                        break;
                    case "fontSizeMultiplier": fontMultiplier = (Integer.parseInt(parameterVal));
                        break;                        
                    case "readcolour_unvaried": UtilityFunctions.getInstance().ReadColour_Unvaried = parameterVal;
                        break;
                    case "readcolour_varied": UtilityFunctions.getInstance().ReadColour_Varied = parameterVal;
                        break;          
                    case "readcolour_insertion": UtilityFunctions.getInstance().ReadColour_Insertion = parameterVal;
                        break;
                    default: System.err.println("Undeclared Parameter "+parameterName);
                }
            }

	} catch (IOException e) {
            e.printStackTrace();
	} 
               
	// Ensure the .scri-bioinf folder exists
	File fldr = new File(filePath, ".scri-bioinf");
	fldr.mkdirs();  

	// Cached reference file
	//ConsensusFileCache.setIndexFile(new File(fldr, "tablet-refs.xml"));
        // This is the file we really want
        File file = new File(fldr, "tablet.xml");
        // So if it exists, just use it
        if (file.exists())
            return file;   
        
        return null;
    }
    
    /*private void CreateTestPanel(){
        
        SetupDebugParameters();
        
        referenceManager[0] = new ReferenceManager(ReferencePath);
        panelManager[0] = new PanelManager();
        fragmentManager = new FragmentManager(DataPath, CachePath);
        String[] refData;
        Map<String,FragmentNode>[] tempFragments;
        int maxReadCount = -1;
        
        //Read Json panel Data
        JSONParser parser = new JSONParser();

	try {
		Object obj = parser.parse(new FileReader(TestFilePath));

		JSONObject jsonObject = (JSONObject) obj;
                tempFragments = new HashMap[jsonObject.size()];
                refData = new String[jsonObject.size()];

                for(Object objHeader:jsonObject.keySet()){
                    
                    //Add to reference
                    String headerVal = (String)objHeader;
                    int index = headerVal.indexOf("-");
                    int colVal = Integer.parseInt(headerVal.substring(0, index));
                    String refVal = headerVal.substring(index+1);
                    refData[colVal] = refVal;
                      
                    //Add to panel
                    tempFragments[colVal] = new HashMap();
                    JSONObject baseJSON = (JSONObject)jsonObject.get(objHeader);
                    int currentReadCount = 0;
                    for(Object baseData:baseJSON.keySet()){
                        String baseVal = (String)baseData;
                        String baseType = baseVal.substring(0, 1);
                        int readCount = Integer.parseInt(baseVal.substring(2)); //@E changed 0-
                        
                        FragmentNode tempNode = new FragmentNode();
                        tempNode.ReadCount = readCount;
                        currentReadCount+=readCount;
                        tempNode.ConnectedFragments = new HashMap();

                        JSONArray connectedJSON = (JSONArray) baseJSON.get(baseData);
                        Iterator<String> iterator = connectedJSON.iterator();

                        while (iterator.hasNext()) {
                            String connectedData = iterator.next();
                            int connectedValIndex = connectedData.indexOf("[");
                            String connectedVal = connectedData.substring(0, connectedValIndex);
                            
                            //Add connected Columns
                            String connectedColumnsRaw = connectedData.substring(connectedValIndex+1,
                                    connectedData.length()-1);
                            String[] connectedCols = connectedColumnsRaw.split(",");     
                            HashSet<Integer> connectedVals = new HashSet();
                            for(String val:connectedCols){
                                connectedVals.add(Integer.parseInt(val));
                            }
                            
                            tempNode.ConnectedFragments.put(connectedVal, connectedVals);
                        }
                        
                        tempFragments[colVal].put(baseType, tempNode); //@E what is temp node
                    }  
                    if(currentReadCount>maxReadCount)
                        maxReadCount = currentReadCount;
                }                
            //Create Panel
            panelManager[0].AddPanel("TestPanel",GridStartX, GridStartY + (PanelSeparation), 
                    FLANK, ColumnWidth, RowHeight,FragmentXOffset,RenderColumns, lineMultiplier, fontMultiplier); 
            panelManager[0].GetPanelFromPhenotype("TestPanel").Fragments = tempFragments;
            panelManager[0].MaxReadCount = maxReadCount;

            //Add to reference Data
            referenceManager[0].AddReference("TestPanel", new ArrayList<String>(Arrays.asList(refData)));
            referenceManager[0].ShiftVals = new HashMap(); //@E for the shift values. Is there not already one?
            referenceManager[0].ShiftVals.put("TestPanel", 0);
            
	} catch (FileNotFoundException e) {
		e.printStackTrace();
	} catch (IOException e) {
		e.printStackTrace();
	} catch (ParseException e) {
		e.printStackTrace();
	}

    }*/
    
    private void SetupDebugParameters(){
        TestFilePath = "C:\\Users\\anwar01\\Documents\\gv15\\test.json";
        WIDTH = 1000;
        HEIGHT = 500;
        GridStartX = 200;
        GridStartY = 100;
        PanelSeparation = 0;
        ColumnWidth = 90;
        RowHeight = 18;
        FragmentXOffset = 20;
        OutputType = "png";
        FLANK = 3;
    }
}