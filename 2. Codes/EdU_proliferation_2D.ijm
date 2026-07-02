// ============================================================
// EdU Proliferation Index — ImageJ/Fiji macro
// Channels: DAPI (w2CSU-405-Em-447) + EdU (w3CSU-640-Em-685-40)
// DAPI segmentation: StarDist | EdU segmentation: "Threshold" mode
// ============================================================

// ---------- Initialization ----------
dir  = getDirectory("Choose the source folder");
list = getFileList(dir);
run("Clear Results");
setBatchMode(true);

outDir  = dir + "Results-Mode" + File.separator;
csvPath = dir + "proliferation_index-Mode.csv";

if (!File.exists(outDir)) File.makeDirectory(outDir);
File.append("Image,Total_DAPI,Positive_EdU,Index_Proliferation_Percent", csvPath);

// ---------- Loop ----------
for (i = 0; i < list.length; i++) {

    name = list[i];
    if (!endsWith(name, "w2CSU-405-Em-447.TIF")) continue;

    prefix   = replace(name, "w2CSU-405-Em-447.TIF", "");
    EdU_name = prefix + "w3CSU-640-Em-685-40.TIF";
    if (!File.exists(dir + EdU_name)) continue;

    print("Analysis : " + prefix);

    // ---------- Open files ----------
    open(dir + name);
    dapi = getTitle();

    open(dir + EdU_name);
    edu = getTitle();


    // ----------  DAPI Segmentation (StarDist) ----------
    selectWindow(dapi);
    run("Subtract Background...", "rolling=50");
    run("Gaussian Blur...", "sigma=2");
    run("Command From Macro",
        "command=[de.csbdresden.stardist.StarDist2D] "
        + "args=['input':'"        + dapi + "', "
        + "'modelChoice':'Versatile (fluorescent nuclei)', "
        + "'normalizeInput':'true', "
        + "'percentileBottom':'1', "
        + "'percentileTop':'99.8', "
        + "'probThresh':'0.5', "
        + "'nmsThresh':'0.4', "
        + "'outputType':'Both', "
        + "'nTiles':'1', "
        + "'excludeBoundary':'2', "
        + "'roiPosition':'Automatic', "
        + "'verbose':'false', "
        + "'showCsbdeepProgress':'false', "
        + "'showProbAndDist':'false'] "
        + ", process=[false]");

    // DAPI binary mask from Label image 
    selectWindow("Label Image");
    setThreshold(1, 65535, "raw");
    run("Convert to Mask");
    saveAs("Tiff", outDir + prefix + "mask_DAPI.TIF");
    dapi_mask = getTitle();
    resetThreshold();
    roiManager("Reset");
    run("Analyze Particles...", "size=100-Infinity display results show=Nothing clear add ");
    nTotal = roiManager("count");
    

	// ---------- Image overlay (Dapi + yellow edges) ----------
    selectWindow(dapi);
    dapi_orig = getTitle();

    run("Grays");
    run("RGB Color");
    setForegroundColor(255, 255, 0);
    roiManager("Show None");
    for (r = 0; r < roiManager("count"); r++) {
        roiManager("Select", r);
        run("Draw");
    }
    run("Select None");
    saveAs("Tiff", outDir + prefix + "overlay_DAPI.TIF");
    close();
    roiManager("Reset");
    
    // ---------- EdU Segmentation (Mode) ----------
    selectWindow(edu);
    run("Set Measurements...", "modal redirect=None decimal=0");
	run("Measure");
	selectWindow("Results");
	modal = getResult("Mode");
	run("Clear Results");
	
    run("Subtract Background...", "rolling=50");
    run("Gaussian Blur...", "sigma=2");
    threshold = 1 * modal;
    print("threshold = "+threshold);
	setThreshold(threshold, 65535);
    setOption("BlackBackground", true);
    run("Convert to Mask");
    run("Fill Holes");
    run("Watershed");
    saveAs("Tiff", outDir + prefix + "mask_EdU.TIF");
    edu_mask = getTitle();
    
    run("Analyze Particles...", "size=100-Infinity display results show=Nothing clear add ");
    nPositive = roiManager("count");

    // ---------- Image overlay (EdU + yellow edges) ----------
    open(dir + EdU_name);
    edu_orig = getTitle();

    run("Grays");
    run("RGB Color");
    setForegroundColor(255, 255, 0);
    roiManager("Show None");
    for (r = 0; r < roiManager("count"); r++) {
        roiManager("Select", r);
        run("Draw");
    }
    run("Select None");
    saveAs("Tiff", outDir + prefix + "overlay_EdU_positive.TIF");
    close();

    // ---------- Calculation and CSV export ----------
    percent  = (nPositive / nTotal) * 100;
    dataLine = prefix + "," + nTotal + "," + nPositive + "," + d2s(percent, 2);
    File.append(dataLine, csvPath);

    // ---------- Cleaning ----------
    roiManager("Reset");
    close(edu_mask);
    close(dapi_mask);
    close(dapi);
    close(edu);
}

// ---------- End ----------
setBatchMode(false);
showMessage("Analyse terminée !\nRésultats dans : " + outDir
            + "\nCSV : " + csvPath);
