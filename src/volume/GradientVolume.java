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

 // This function computes the magnitude of the local gradient vector at each data point
    //and stores it into the data attribute. 
    private void compute() {

        for (int i=0; i<data.length; i++) {
            data[i] = zero;
        }
        
        //computes the gradient in each axes seperately
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
    }
    	
    //This function linearly interpolates gradient vector g0 and g1 given the factor (t).
    private VoxelGradient interpolate(VoxelGradient g0, VoxelGradient g1, float factor) {

        
        float x = g1.x*(1 - factor) + factor*g0.x;
        float y = g1.y*(1 - factor) + factor*g0.y;
        float z = g1.z*(1 - factor) + factor*g0.z;
        
        //result.x = g1.x*(1-factor) + g0.x*factor;
        //result.y = g1.y*(1-factor) + g0.y*factor;
        //result.z = g1.z*(1-factor) + g0.z*factor;
        //result.mag = (float) Math.sqrt(result.x * result.x + result.y * result.y + result.z * result.z);

        VoxelGradient result = new VoxelGradient(x,y,z);
        return result;
    }

    //We first interpolate between the lower x axis (x0, x2), then the upper x axis (x1, x3), and then
    //we interpolate these two results.
    /**
     * 
     * @param g0 - value at "left" vertex on lower x axis
     * @param g1 - value at "left" vertex on upper x axis
     * @param g2 - value at "right" vertex on lower x axis
     * @param g3 - value at "left" vertex on upper x axis
     * @param factor_1 - factor between the 2 x axes
     * @param factor_2 - factor between the interpolated values
     * @return 
     */
 
    private VoxelGradient interpolate2D(VoxelGradient g0, VoxelGradient g1, VoxelGradient g2, VoxelGradient g3, float factor_1, float factor_2){
        VoxelGradient ans1 = interpolate(g0, g2, factor_1);
        VoxelGradient ans2 = interpolate(g1, g3, factor_1);
        VoxelGradient result =  interpolate(ans1, ans2, factor_2);
        
        return result;
    }
   /**
    	//This function linearly interpolates gradient vector g0 and g1 given the factor (t) 
    //the resut is given at result. You can use it to tri-linearly interpolate the gradient 
	private void interpolate(VoxelGradient g0, VoxelGradient g1, float factor, VoxelGradient result) {
        result.x = (1.0f - factor)*g0.x + factor*g1.x;
        result.y = (1.0f - factor)*g0.y + factor*g1.y;
        result.z = (1.0f - factor)*g0.z + factor*g1.z;
        result.mag = (float) Math.sqrt(result.x * result.x + result.y * result.y + result.z * result.z);
    }**/
        
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
        
        //get the gradients at each vertices 
        VoxelGradient c000 = getGradient(x0,y0,z0);
        VoxelGradient c001 = getGradient(x0,y0,z1);
        VoxelGradient c010 = getGradient(x0,y1,z0);
        VoxelGradient c100 = getGradient(x1,y0,z0);
        VoxelGradient c011 = getGradient(x0,y1,z1);
        VoxelGradient c101 = getGradient(x1,y0,z1);
        VoxelGradient c110 = getGradient(x1,y1,z0);
        VoxelGradient c111 = getGradient(x1,y1,z1);
        
        VoxelGradient ans1 = interpolate2D(c000, c010, c100, c110, (float) (x - x0), (float) (y - y0));
        VoxelGradient ans2 = interpolate2D(c001, c011, c101, c111, (float) (x - x0), (float) (y - y0));
        VoxelGradient result = interpolate(ans1, ans2, (float) z - z0);
        
        //compute the interpolated gradients
        
        return result;
    }
    /**
        public VoxelGradient getGradient(double[] coord) {

        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        float fac_x = (float) coord[0] - x;
        float fac_y = (float) coord[1] - y;
        float fac_z = (float) coord[2] - z;

        VoxelGradient t0 = new VoxelGradient();
        interpolate(getGradient(x, y, z), getGradient(x+1, y, z), fac_x, t0);
        VoxelGradient t1 = new VoxelGradient();
        interpolate(getGradient(x, y+1, z), getGradient(x+1, y+1, z), fac_x, t1);
        VoxelGradient t2 = new VoxelGradient();
        interpolate(getGradient(x, y, z+1), getGradient(x+1, y, z+1), fac_x, t2);
        VoxelGradient t3 = new VoxelGradient();
        interpolate(getGradient(x, y+1, z+1), getGradient(x+1, y+1, z+1), fac_x, t3);

        VoxelGradient t4 = new VoxelGradient();
        interpolate(t0, t1, fac_y, t4);
        VoxelGradient t5 = new VoxelGradient();
        interpolate(t2, t3, fac_y, t5);
        
        VoxelGradient t6 = new VoxelGradient();
        interpolate(t4, t5, fac_z, t6);
        return t6;

    }**/
    
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
