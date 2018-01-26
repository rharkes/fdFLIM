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
package Classes;
        
public class PhiMod {
    private final double[] sinvec;
    private final double[] cosvec;
    private final int nr_phases;
    public double modulation;
    public double phase;
    public double GoodFit;
    public double mean;
    
    
    public PhiMod(int phases) {
        nr_phases = phases;
        sinvec = new double[phases];
        cosvec = new double[phases];
        for (int i=0;i<nr_phases;i++){
            double frac = (double)i / (double)phases;
            sinvec[i] = Math.sin(2*Math.PI*frac);
            cosvec[i] = Math.cos(2*Math.PI*frac);
        }
    }
    public void calc(double data[]){
        if (data.length!=nr_phases){throw new java.lang.Error("Wrong nr phases");}
        double F0 = 0;
        double Fs = 0;
        double Fc = 0;
        for (int i=0;i<nr_phases;i++){
            F0 += data[i];
            Fs += (sinvec[i]*data[i]);
            Fc += (cosvec[i]*data[i]);
        }
        mean = F0/nr_phases;
        Fs = Fs/nr_phases;
        Fc = Fc/nr_phases;
        modulation = 2*Math.sqrt((Fs*Fs)+(Fc*Fc))/mean;
        phase = Math.atan2(Fc, Fs);
        GoodFit=1;
    }
}
