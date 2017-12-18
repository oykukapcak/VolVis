/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class GradientVolume {

	
	
	//////////////////////////////////////////////////////////////////////
	///////////////// TO BE IMPLEMENTED //////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
	//Compute the gradient of contained in the volume attribute and save it into the data attribute
	//
	//This is a lengthy computation and is performed only once (have a look at the constructor GradientVolume) 

    private void compute() {

        for (int i=0; i<data.length; i++) {
            data[i] = zero;
        }
        
        for (int z=1; z<dimZ-1; z++) {
            for (int y=1; y<dimY-1; y++) {
                for (int x=1; x<dimX-1; x++) {
                    float gx = (volume.getVoxel(x+1, y, z) - volume.getVoxel(x-1, y, z))/2.0f;
                    float gy = (volume.getVoxel(x, y+1, z) - volume.getVoxel(x, y-1, z))/2.0f;
                    float gz = (volume.getVoxel(x, y, z+1) - volume.getVoxel(x, y, z-1))/2.0f;
                    VoxelGradient grad = new VoxelGradient(gx, gy, gz);
                    setGradient(x, y, z, grad);
                }
            }
        }
        
    }
    	
	//This function linearly interpolates gradient vector g0 and g1 given the factor (t) 
    //the resut is given at result. You can use it to tri-linearly interpolate the gradient 
	private void interpolate(VoxelGradient g0, VoxelGradient g1, float factor, VoxelGradient result) {
        result.x = (1.0f - factor)*g0.x + factor*g1.x;
        result.y = (1.0f - factor)*g0.y + factor*g1.y;
        result.z = (1.0f - factor)*g0.z + factor*g1.z;
        result.mag = (float) Math.sqrt(result.x * result.x + result.y * result.y + result.z * result.z);
    }
	
	//Do NOT modify this function
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

    }
    
    //Returns the maximum gradient magnitude
    //
    //The data array contains all the gradients, in this function you have to return the maximum magnitude of the vectors in data[] 
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
    	
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
	
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

  
	//Do NOT modify this function
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
