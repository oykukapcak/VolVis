/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 *  Modified by Anna Vilanova
 */
public class GradientVolume {

	
	
	//////////////////////////////////////////////////////////////////////
	///////////////// TO BE IMPLEMENTED //////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
	//Compute the gradient of contained in the volume attribute and save it into the data attribute
	//
	//This is a lengthy computation and is performed only once (have a look at the constructor GradientVolume) 

    // You need to implement this function
    private void compute() {

        for (int i=0; i<data.length; i++) {
            data[i] = zero;
        }
        for(int x = 1; x<dimX-1; x++){
            for(int y = 1; y<dimY-1; y++){
                for(int z = 1; z<dimZ-1; z++){
                    short x1 = volume.getVoxel(x-1,y,z);
                    short x2 = volume.getVoxel(x+1,y,z);
                    short y1 = volume.getVoxel(x,y-1,z);
                    short y2 = volume.getVoxel(x,y+1,z);
                    short z1 = volume.getVoxel(x,y,z-1);
                    short z2 = volume.getVoxel(x,y,z+1);
                    
                    setGradient(x,y,z, new VoxelGradient((float)(x2-x1)/2, (float)(y2-y1)/2, (float)(z2-z1)/2));
                }
            }
        }
        // to be implemented
    }
    	
    //You need to implement this function
    //This function linearly interpolates gradient vector g0 and g1 given the factor (t) 
    //the resut is given at result. You can use it to tri-linearly interpolate the gradient
    private void interpolate(VoxelGradient g0, VoxelGradient g1, float factor, VoxelGradient result) {

        float x = g1.x*(1 - factor)/factor + g0.x;
        float y = g1.y*(1 - factor)/factor + g0.y;
        float z = g1.z*(1 - factor)/factor + g0.z;
  
        result = new VoxelGradient(x,y,z); 
    }

    //g0 - value at 0,0, g1 - value at 0,1, g2 - value at 1,0, g3 - value at 1,1
    private void interpolate2D(VoxelGradient g0, VoxelGradient g1, VoxelGradient g2, VoxelGradient g3, float factor_1, float factor_2, VoxelGradient result){
        VoxelGradient ans1 = null;
        VoxelGradient ans2 = null;
        VoxelGradient ans3 = null;
        interpolate(g0, g2, factor_1, ans1);
        interpolate(g1, g3, factor_1, ans2);
        interpolate(ans1, ans2, factor_2, ans3);
        result = ans3;
    }
    
    // You need to implement this function
    // This function returns the gradient at position coord using trilinear interpolation
    public VoxelGradient getGradient(double[] coord) {

        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x0 = (int) Math.floor(coord[0]); 
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;
        
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        
        VoxelGradient c000 = getGradient(x0,y0,z0);
        VoxelGradient c001 = getGradient(x0,y0,z1);
        VoxelGradient c010 = getGradient(x0,y1,z0);
        VoxelGradient c100 = getGradient(x1,y0,z0);
        VoxelGradient c011 = getGradient(x0,y1,z1);
        VoxelGradient c101 = getGradient(x1,y0,z1);
        VoxelGradient c110 = getGradient(x1,y1,z0);
        VoxelGradient c111 = getGradient(x1,y1,z1);
        
        VoxelGradient ans1 = null;
        VoxelGradient ans2 = null;
        VoxelGradient ans3 = null;
        
        interpolate2D(c000, c010, c100, c110, (float) (x - x0), (float) (y - y0), ans1);
        interpolate2D(c001, c011, c101, c111, (float) (x - x0), (float) (y - y0), ans2);
        interpolate(ans1, ans2, (float) z - z0, ans3);
        return ans3;
    }
    
    
    
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
	
    //Returns the maximum gradient magnitude
    //
    //The data array contains all the gradients, in this function you have to return the maximum magnitude of the vectors in data[] 
 
    //Do NOT modify this function
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    	
	
	
	//Do NOT modify this function
	public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

	//Do NOT modify this function
	public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

  
	//Do NOT modify this function: Basically calculates the Nearest Neighbor interpolation for the gradient
    public VoxelGradient getGradient2(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
        return getGradient(x, y, z);
    }

	//Do NOT modify this function
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

	//Do NOT modify this function
    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }
    
	//Do NOT modify this function
    public VoxelGradient getVoxel(int i) {
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

	//Do NOT modify this attributes
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
    
    //If needed add new attributes here:
}
