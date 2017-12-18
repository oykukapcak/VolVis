/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author michel
 */
public class Volume {
    
	//////////////////////////////////////////////////////////////////////
	///////////////// TO BE IMPLEMENTED //////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
	//This function linearly interpolates the value g0 and g1 given the factor (t) 
    //the result is returned. You can use it to tri-linearly interpolate the values 
	private float interpolate(float g0, float g1, float factor) {
        float result = (1.0f - factor)*g0 + factor*g1;
        return result; 
    }
	
	//You have to implement the trilinear interpolation of the volume
	//First implement the interpolated funzion above
	public short getVoxelInterpolate(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.floor(coord[0]); 
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        float fac_x = (float) coord[0] - x;
        float fac_y = (float) coord[1] - y;
        float fac_z = (float) coord[2] - z;

        float t0 = interpolate(getVoxel(x, y, z), getVoxel(x+1, y, z), fac_x);
        float t1 = interpolate(getVoxel(x, y+1, z), getVoxel(x+1, y+1, z), fac_x);
        float t2 = interpolate(getVoxel(x, y, z+1), getVoxel(x+1, y, z+1), fac_x);
        float t3 = interpolate(getVoxel(x, y+1, z+1), getVoxel(x+1, y+1, z+1), fac_x);
        float t4 = interpolate(t0, t1, fac_y);
        float t5 = interpolate(t2, t3, fac_y);
        float t6 = interpolate(t4, t5, fac_z);
        
        return (short)t6; 
    }
		
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	//Do NOT modify this function
	public short getVoxelNN(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-1) || coord[1] < 0 || coord[1] > (dimY-1)
                || coord[2] < 0 || coord[2] > (dimZ-1)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.round(coord[0]); 
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
    
        return getVoxel(x, y, z);
    }
	
	//Do NOT modify this function
    public Volume(int xd, int yd, int zd) {
        data = new short[xd*yd*zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }
	//Do NOT modify this function
    public Volume(File file) {
        
        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }
        
    }
    
	//Do NOT modify this function
    public short getVoxel(int x, int y, int z) {
    	int i = x + dimX*(y + dimY * z);
        return data[i];
    }
    
	//Do NOT modify this function
    public void setVoxel(int x, int y, int z, short value) {
    	int i = x + dimX*(y + dimY * z);
        data[i] = value;
    }
    
	//Do NOT modify this function
    public void setVoxel(int i, short value) {
        data[i] = value;
    }
    
	//Do NOT modify this function
    public short getVoxel(int i) {
        return data[i];
    }
    
	//Do NOT modify this function
    public int getDimX() {
        return dimX;
    }
    
	//Do NOT modify this function
    public int getDimY() {
        return dimY;
    }
    
	//Do NOT modify this function
    public int getDimZ() {
        return dimZ;
    }

	//Do NOT modify this function
    public short getMinimum() {
        short minimum = data[0];
        for (int i=0; i<data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }
    
	//Do NOT modify this function
    public short getMaximum() {
        short maximum = data[0];
        for (int i=0; i<data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }
 
	//Do NOT modify this function
    public int[] getHistogram() {
        return histogram;
    }
    
	//Do NOT modify this function
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i=0; i<data.length; i++) {
            histogram[data[i]]++;
        }
    }
    
	//Do NOT modify these attributes
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;
}
