/* 
 * Copyright (C) 2018 Rolf Harkes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package Main;

import Classes.*;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import net.imagej.ImageJ;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/**
 * Calculates lifetimes from .fli files.
 */
@Plugin(type = Command.class, headless = true,
        menuPath = "Plugins>fdFLIM")
public class fdFLIM implements Command, Previewable {

    @Parameter
    private LogService log;

    @Parameter
    private StatusService statusService;
    
    @Parameter(label = "Select Sample Image", description = "Sample Image")
    private ImagePlus Image1;
    
    @Parameter(label = "Only Phi and Mod", description = "Calculate only phase shift and modulation")
    private Boolean BoolPhiMod;
    
    @Parameter(label = "Goodness of Fit", description = "Calculate the Goodness of Fit parameter")
    private Boolean GF;
    
    @Parameter(label = "Select Reference Image", description = "Reference Image")
    private ImagePlus Image2;
    
    @Parameter(label = "Reference lifetime", description = "Reference lifetime")
    private Double Tau_ref;
    
    @Parameter(label = "Measurement Frequency", description = "Measurement Frequency")
    private Double Freq;

    public static void main(final String... args) throws Exception {
        // Launch ImageJ as usual.
        final ImageJ ij = net.imagej.Main.launch(args);
        // Launch the command.
        ij.command().run(fdFLIM.class, true);
    }

    @Override
    public void run() {
        ImageStack stack1 = Image1.getStack();
        ImageStack result1 = getPhiMod(stack1);
        if (BoolPhiMod){
            ImagePlus ResIm = new ImagePlus("Phase Modulation Mean",result1);
            ResIm.show();
            return;
        }
        ImageStack stack2 = Image2.getStack();
        ImageStack result2 = getPhiMod(stack2);
        //get lifetimes
        if (Freq<1E3){Freq=Freq*1E6;}
        if (Tau_ref>1E-5){Tau_ref=Tau_ref*1E-9;}
        double Omega = 2*Math.PI*Freq;
        double phi_ref = Math.atan(Omega*Tau_ref);
        double mod_ref = Math.sqrt(1/(Math.pow(Omega*Tau_ref,2)+1));
        ImageStack lifetimes = getLifetime(result2,result1,phi_ref,mod_ref,Omega);
        ImagePlus ResIm = new ImagePlus("Lifetimes",lifetimes);
        ResIm.show();
    }

    @Override
    public void cancel() {
        log.info("Cancelled");
    }

    @Override
    public void preview() {
        log.info("Preview");
    }
    private ImageStack getLifetime(ImageStack Ref, ImageStack Sam, double phi_ref,double mod_ref,double Omega){
        int h = Ref.getHeight();
        int w = Ref.getWidth();
        ImageStack result = new ImageStack(w,h);
        result.addSlice("Lifetime from Phase",new ByteProcessor(w,h));
        result.addSlice("Lifetime from Modulation",new ByteProcessor(w,h));
        result.addSlice("Mean",new ByteProcessor(w,h));
        if (GF){result.addSlice("Goodness of Fit",new ByteProcessor(w,h));}
        result = result.convertToFloat();
        for (int x = 0;x<w;x++){
            for (int y=0;y<h;y++){
                //phase
                double phi = Sam.getVoxel(x,y,0) - (Ref.getVoxel(x, y, 0)-phi_ref);
                double mod = (Sam.getVoxel(x,y,1)*mod_ref)/Ref.getVoxel(x, y, 1);
                result.setVoxel(x, y, 0, 1E9*Math.tan(phi)/Omega);
                result.setVoxel(x, y, 1, 1E9*Math.sqrt((1/(mod*mod))-1)/Omega);
                result.setVoxel(x, y, 2, Sam.getVoxel(x,y,2));
                if (GF){
                    result.setVoxel(x, y, 3, Math.sqrt(Math.pow(Sam.getVoxel(x,y,4),2)+Math.pow(Ref.getVoxel(x,y,4),2)));
                }
            }
        }
        return result;
    }
    private ImageStack getPhiMod(ImageStack input){
        int phases = input.getSize();
        int h = input.getHeight();
        int w = input.getWidth();
        PhiMod PM = new PhiMod(phases);
        PM.doGF = GF;
        ImageStack result = new ImageStack(w,h);
        result.addSlice("Phase",new ByteProcessor(w,h));
        result.addSlice("Modulation",new ByteProcessor(w,h));
        result.addSlice("Mean",new ByteProcessor(w,h));
        if (PM.doGF){result.addSlice("Goodness of Fit",new ByteProcessor(w,h));}
        result = result.convertToFloat();
        double[] data = new double[phases];
        for (int x = 0;x<w;x++){
            for (int y=0;y<h;y++){
                for (int ph=0;ph<phases;ph++){
                    data[ph] = input.getVoxel(x, y, ph);
                }
                PM.calc(data);
                result.setVoxel(x, y, 0, PM.phase);
                result.setVoxel(x, y, 1, PM.modulation);
                result.setVoxel(x, y, 2, PM.mean);
                if (PM.doGF){result.setVoxel(x, y, 2, PM.GoodFit);}
            }
        }
        return result;
    }
}