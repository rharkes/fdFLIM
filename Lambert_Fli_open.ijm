@File(label = "Sample File", style = "file") Sample
@Float(label = "Frequency", style = "spinner", value=40) freq
@File(label = "Reference File", style = "file") Reference
@Float(label = "Reference lifetime", style = "spinner", value=5.4) tau

//Macro to open a reference and sample file and calculate lifetimes from phase and modulation.

run("Close All");
OpenFile(Reference)
ref = getTitle();
OpenFile(Sample)
sam = getTitle();
run("fdFLIM", "image1=["+sam+"] boolphimod=false image2=["+ref+"] tau_ref="+tau+" freq="+freq);
run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames=1 display=Color");
close(sam);close(ref);
Stack.setChannel(3);
setAutoThreshold("Otsu dark");
run("Create Mask");
run("Divide...", "value=255");
run("32-bit");
changeValues(0,0,NaN);
imageCalculator("Multiply 32-bit stack", "Lifetimes","mask");
close("mask");
Stack.setChannel(1);run("Rainbow RGB");setMinAndMax(1, 4);
Stack.setChannel(2);run("Rainbow RGB");setMinAndMax(1, 4);
Stack.setChannel(3);run("Grays");run("Enhance Contrast", "saturated=0.35");


function OpenFile(file) {
	run("Bio-Formats", "open=["+file+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
	image = getTitle();
	run("Bio-Formats", "open=["+file+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_2");
	background = getTitle();
	imageCalculator("Subtract create 32-bit stack",image,background);
	close(image)
	close(background)
}
