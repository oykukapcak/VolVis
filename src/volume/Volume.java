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
 *  Modified by Anna Vilanova
 */
public class Volume {
    
	//////////////////////////////////////////////////////////////////////
	///////////////// TO BE IMPLEMENTED //////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
    //This function linearly interpolates the value x0 and x1 given the factor (t) 
    //the result is returned. You can use it to tri-linearly interpolate the values 
	private float interpolate1D(float g0, float g1, float factor) {
        float result=0;
        result = factor*(g1 - g0)+ g0;
        return result; 
    }
       //We first interpolate between the lower x axis (x0, x2), then the upper x axis (x1, x3), and then
       //we interpolate by
        /**
         * 
         * @param x0 - "left" vertex on lower x axis
         * @param x1 - "left" vertex on upper x axis
         * @param x2 - "right" vertex on lower x axis
         * @param x3 - "left" vertex on upper x axis
         * @param factor_x - x-point between the 2 x axes
         * @param factor_y - y-point between the interpolated values
         * @return 
         */
        private float interpolate2D(float x0, float x1, float x2, float x3, float factor_x, float factor_y){
            float x_result1 = interpolate1D(x0, x2, factor_x);
            float x_result2 = interpolate1D(x1, x3, factor_x);
            float result = interpolate1D(x_result1, x_result2, factor_y);
            return result;
        }
	
	//You have to implement the trilinear interpolation of the volume
	//First implement the interpolated function above
        // At the moment the function does takes just the lowest voxel value
        // to trilinear interpolation
	public short getVoxelLinearInterpolate(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x0 = (int) Math.floor(coord[0]); 
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;
        
        double x = coord[0] - x0;
        double y = coord[1] - y0;
        double z = coord[2] - z0;
        
        // get the vertices of the cube
        short c000 = getVoxel(x0,y0,z0);
        short c001 = getVoxel(x0,y0,z1);
        short c010 = getVoxel(x0,y1,z0);
        short c100 = getVoxel(x1,y0,z0);
        short c011 = getVoxel(x0,y1,z1);
        short c101 = getVoxel(x1,y0,z1);
        short c110 = getVoxel(x1,y1,z0);
        short c111 = getVoxel(x1,y1,z1);
        
        float result_z0 = interpolate2D(c000, c010, c100, c110, (float)x, (float)y);
        float result_z1 = interpolate2D(c001, c011, c101, c111, (float)x, (float)y);
        float result = interpolate1D(result_z0, result_z1, (float)z);
            
        return (short)result;
    }
		
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	//Do NOT modify this function
        // This function is an example and does a nearest neighbour interpolation
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
