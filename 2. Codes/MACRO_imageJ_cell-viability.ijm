// Conditions of the experiments to FILL
model = "BC152";
tumor = 4;
day = 0;
concentration = 0; // in µM
chip = "A";

// Set replica based on chip
if (chip == "A") {
    replica = 1;
} else if (chip == "B") {
    replica = 2;
} else if (chip == "C") {
    replica = 3;
}

// Set drug based on concentration
if (concentration == 0) {
    drug = "Control";
} else {
    drug = "Carboplatin";
}

// Choice of the repertory for saving the .txt file
output_path = getDirectory("Choose where to save results file:") + "results_" + model + "_" + drug + "_Chip" + chip + "_Day" + day + "_Conc" + concentration + ".txt";

// Création du fichier et écriture de l'en-tête
file = File.open(output_path);
File.append("Drug\tTumor\tReplica\tDay\tConcentration\tN_Dead\tN_Live", output_path);


// Ask the user for the directory containing the images to process
directory = getDirectory("Select the directory containing the images:");

// Get a list of files in the directory
fileList = getFileList(directory);

// Process each image in the directory
for (j = 0; j < fileList.length/4; j++) { //fileList.length/4
	// Get the file path of the current image
	filePath_BF = directory + fileList[4*j+1];
	filePath_L = directory + fileList[4*j+2];
	filePath_D = directory + fileList[4*j+3];
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Process dead image, counting of dead cells
	open(filePath_D);
	rename("Dead");
	setOption("ScaleConversions", true);
	//run("8-bit");

	run("Set Measurements...", "modal redirect=None decimal=0");
	run("Measure");
	selectWindow("Results");
	modal = getResult("Mode");	
	
	run("Gaussian Blur...", "sigma=2");	
	
		//Binarization
	selectWindow("Dead");
	
	threshold = 3 * modal;
	setThreshold(threshold, 65535);
	
	//setOption("BlackBackground", false);
	run("Convert to Mask");
		//Binarization post-processing
	run("Watershed");
		//Remove small objects
	run("Analyze Particles...", "size=0-Infinity show=Masks");
	selectWindow("Mask of Dead");
		//Count Number of Dead cells
	run("Set Measurements...", "area redirect=None decimal=0");
	run("Analyze Particles...", "size=30-Infinity display results show=Overlay display");
	N_dead = nResults;
	
	print("Experience" + j+1);
	fluorescentArea = 0;
for (i = 0; i < nResults; i++) {
    fluorescentArea += getResult("Area", i);
}

//"Dead fluorescence area
print(fluorescentArea);
	selectWindow("Results");
	run("Close");
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Process Live image, countin of live cells
	open(filePath_L);
	rename("Live");
	setOption("ScaleConversions", true);
	//run("8-bit");

	run("Set Measurements...", "modal redirect=None decimal=0");
	run("Measure");
	selectWindow("Results");
	modal = getResult("Mode");	
		
	run("Gaussian Blur...", "sigma=2");		
		
		//Binarization
	selectWindow("Live");
	
	threshold = 3 * modal;
	setThreshold(threshold, 65535);
	
	//setOption("BlackBackground", false);
	run("Convert to Mask");
		//Binrization post-processing
	run("Watershed");
		//Remove small objects
	run("Analyze Particles...", "size=0-Infinity show=Masks");
	selectWindow("Mask of Live");
		//Count Number of Dead cells
	run("Set Measurements...", "area redirect=None decimal=0");
	run("Analyze Particles...", "size=65-Infinity display results show=Overlay display");
	N_Live = nResults;
	
	fluorescentArea = 0;
for (i = 0; i < nResults; i++) {
    fluorescentArea += getResult("Area", i);
}
//Live fluorescence area:
print(fluorescentArea);
	selectWindow("Results");
	run("Close");
	
	// Print the numbers and close the current windows
	print(N_dead);	
	print(N_Live);
	selectWindow("Dead");
	run("Close");
	selectWindow("Live");
	run("Close");
	selectWindow("Mask of Live");
	run("Close");
	selectWindow("Mask of Dead");
	run("Close");
	
	// Call in the loop the experimental conditions displayed outside of the loop
	print("drug = " + drug);
	print("tumor = " + tumor);
	print("replica = " + replica);
	print("day = " + day);
	print("concentration = " + concentration);
	print("N_Live = " + N_Live);
	print("N_dead = " + N_dead);
	
	// Define the order of the line and create the .txt file
	line = "" + drug + "\t" + tumor + "\t" + replica + "\t" + day + "\t" + concentration + "\t" + N_dead + "\t" + N_Live;
	print("Écriture : " + line);
	File.append(line, output_path);
		
}