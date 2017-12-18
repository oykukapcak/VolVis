/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

	
	//////////////////////////////////////////////////////////////////////
	///////////////// TO BE IMPLEMENTED //////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	//in this function we update the "image" attribute using the slicing technique
    void slicer(double[] viewMatrix) {
	    // we start by clearing the image
	    resetImage();
	
	    // vector uVec and vVec define the view plane, 
	    // perpendicular to the view vector viewVec which is going from the view point towards the object
	    double[] viewVec = new double[3];
	    double[] uVec = new double[3];
	    double[] vVec = new double[3];
	    getPlane(viewVec,uVec,vVec);
	  
	    // compute the volume center
	    double[] volumeCenter = new double[3];
	    computeVolumeCenter(volumeCenter);
	    
	    // Here will be stored the 3D coordinates of every pixel in the plane 
	    double[] pixelCoord = new double[3];
	   	
	    // sample on a plane through the origin of the volume data
	    double max = volume.getMaximum();
	    TFColor voxelColor = new TFColor();
	    //Iterate on every pixel
	    for (int j = 0; j < image.getHeight(); j++) {
	        for (int i = 0; i < image.getWidth(); i++) {
	        	//pixelCoord now contains the 3D coordinates for pixel (i,j)
	        	computePixelCoordinates(pixelCoord,volumeCenter,uVec,vVec,i,j);
		            
	            //pixelCoord now contains the 3D coordinates of the pixels
	            //we now have to get the value for the in the 3D volume for the pixel
	            //we can use a nearest neighbor implementation like this:
	            int val = volume.getVoxelNN(pixelCoord);

	            		
	            //you have to implement the function getVoxelInterpolated in Volume.java
	            //in order to complete the assignment
	            //int val = volume.getVoxelInterpolate(pixelCoord); //and then use this line
	            
	            
	            // Map the intensity to a grey value by linear scaling
	            voxelColor.r = val/max;
	            voxelColor.g = voxelColor.r;
	            voxelColor.b = voxelColor.r;
	            
	            // the following linke makes intensity 0 completely transparent and the rest opaque
	            // Alternatively, apply the transfer function to obtain a color using the tFunc attribute
	            // voxelColor = tFunc.getColor(val);
	            voxelColor.a = val > 0 ? 1.0 : 0.0;  
	            
	            
	            //BufferedImage expects a pixel color packed as ARGB in an int
	            //use the function computeImageColor to convert your double color in the range 0-1 to the format need by the image
	            int pixelColor = computeImageColor(voxelColor.r,voxelColor.g,voxelColor.b,voxelColor.a);
	            image.setRGB(i, j, pixelColor);
	        }
	    }
	}
    

    //Implementation of MIP
    //you should return the color assigned to a ray given it's starting point and the direction.
    //ray must be sampled "nrSamples" times, with a sample step of sampleStep
    int traceRayMIP(double[] startPoint, double[] direction, double sampleStep, int nrSamples) {

    	//Temporary vector used for the sampling
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, startPoint[0], startPoint[1], startPoint[2]);
        
        double accumulator = 0;
        for (int k = 0; k < nrSamples; k++) {
            currentPos[0] = startPoint[0] + k * sampleStep * direction[0];
            currentPos[1] = startPoint[1] + k * sampleStep * direction[1];
            currentPos[2] = startPoint[2] + k * sampleStep * direction[2];
            
            double value = volume.getVoxelInterpolate(currentPos)/255.;
            if (value > accumulator) {
                accumulator = value;
            }
        }

        double alpha;
        double r, g, b;
        if (accumulator > 0.0) { // if the maximum = 0 make the voxel transparent
            alpha = 1.0;
        } else {
            alpha = 0.0;
        }
        r = g = b = accumulator;
        int color = computeImageColor(r,g,b,alpha);
        return color;
    }
    
    //Implementation of a ray tracer given the entry and exit points
    //can you adapt the code of the MIP in order to compute your own nrSamples?
    int traceRay(double[] entryPoint, double[] exitPoint, double[] viewVec, double sampleStep) {
    	//compute the increment and the number of samples
        double[] increments = new double[3];
        VectorMath.setVector(increments, -viewVec[0] * sampleStep, -viewVec[1] * sampleStep, -viewVec[2] * sampleStep);
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);

        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, entryPoint[0], entryPoint[1], entryPoint[2]);
        
        double accumulator = 0;
        do {
            double value = volume.getVoxelInterpolate(currentPos)/255.;
            if (value > accumulator) {
                accumulator = value;
            }
            for (int i = 0; i < 3; i++) {
                currentPos[i] += increments[i];
            }
            nrSamples--;
        } while (nrSamples > 0);

        double alpha;
        double r, g, b;
        if (accumulator > 0.0) { // if the maximum = 0 make the voxel transparent
            alpha = 1.0;
        } else {
            alpha = 0.0;
        }
        r = g = b = accumulator;
        int color = computeImageColor(r,g,b,alpha);
        return color;
    }
    
    
    
    int traceRayComposite(double[] entryPoint, double[] exitPoint, double[] viewVec, double sampleStep) {
        double[] lightVector = new double[3];
        double[] halfVector = new double[3];
        //the light vector is directed toward the view point (which is the source of the light)
        //half vector is used to speed up the phong shading computation
        getLightVector(lightVector,halfVector,viewVec);
        
        //compute increments along the ray for each step
        double[] increments = new double[3];
        computeIncrements(increments, viewVec, sampleStep);
        //compute the number of steps
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);

        //I use a front-to-back composition, therefore the current position is initialized as the exit point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, exitPoint[0], exitPoint[1], exitPoint[2]);

        //Initialization of the colors as floating point values
        double r, g, b;
        r = g = b = 0.0;
        double alpha = 0.0;
        double opacity = 0;
        
        TFColor voxel_color = new TFColor();
        do {
        	//Gets the value and the gradient in the current position
            int value = volume.getVoxelInterpolate(currentPos);
            VoxelGradient gradient = gradients.getGradient(currentPos);

            //Updates the color and the opacity based on the current selection
            if (compositingMode) {
                voxel_color = tFunc.getColor(value);
                opacity = voxel_color.a;    
            }
            if (tf2dMode) {
                voxel_color = tfEditor2D.triangleWidget.color;
                opacity = tfEditor2D.triangleWidget.color.a;            
                opacity *= computeLevoyOpacity(tfEditor2D.triangleWidget.baseIntensity, 
                    tfEditor2D.triangleWidget.radius, value, gradient.mag);
            }
            if (shadingMode) {
                if (opacity > 0.0) {
                    voxel_color = computePhongShading(voxel_color, gradient, lightVector, halfVector);
                }
            }
            
            // Compute the composition with the back-to-front algorithm
            r = opacity * voxel_color.r + (1.0 - opacity) * r;
            g = opacity * voxel_color.g + (1.0 - opacity) * g;
            b = opacity * voxel_color.b + (1.0 - opacity) * b;
            alpha = opacity + (1.0 - opacity) * alpha;
            


            // front-to-back; note: change sign of increments and entry/exit
           /*
             r += voxel_color.a * voxel_color.r * (1.0 - alpha);
             g += voxel_color.a * voxel_color.g * (1.0 - alpha);
             b += voxel_color.a * voxel_color.b * (1.0 - alpha);
             alpha += (1.0-alpha)*voxel_color.a;
             */

            //update the current position
            for (int i = 0; i < 3; i++) {
                currentPos[i] += increments[i];
            }
            nrSamples--;
        } while (nrSamples > 0);
        
        //computes the color
        int color = computeImageColor(r,g,b,alpha);
        return color;
    }
    
    void raycast(double[] viewMatrix) {

    	//data allocation
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];

        // ray parameters
        int increment;
        double sampleStep;
        if (interactiveMode) {
            increment = 2;
            sampleStep = 4.0;
        } else {
            increment = 1;
            sampleStep = 1.0;
        }
        
        // reset the image to black
        resetImage();
        // compute the view plane and the view vector that is used to compute the entry and exit point of the ray
        getPlane(viewVec,uVec,vVec);

        // ray computation for each pixel
        for (int j = 0; j < image.getHeight(); j += increment) {
            for (int i = 0; i < image.getWidth(); i += increment) {
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
            	computePixelCoordinatesBehind(pixelCoord,viewVec,uVec,vVec,i,j);
            	// compute the entry and exit point of the ray
                computeEntryAndExit(pixelCoord, viewVec, entryPoint, exitPoint);
                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                    int val = 0;
                    if (compositingMode || tf2dMode) {
                        val = traceRayComposite(entryPoint, exitPoint, viewVec, sampleStep);
                    } else if (mipMode) {
                        val = traceRay(entryPoint, exitPoint, viewVec, sampleStep);
                    }
                    for (int ii = i; ii < i + increment; ii++) {
                        for (int jj = j; jj < j + increment; jj++) {
                            image.setRGB(ii, jj, val);
                        }
                    }
                }

            }
        }
    }
		
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    private boolean mipMode = false;
    private boolean slicerMode = true;
    private boolean compositingMode = false;
    private boolean tf2dMode = false;
    private boolean shadingMode = false;

    //Do NOT modify this function
    int computeImageColor(double r, double g, double b, double a){
		int c_alpha = 	a <= 1.0 ? (int) Math.floor(a * 255) : 255;
        int c_red = 	r <= 1.0 ? (int) Math.floor(r * 255) : 255;
        int c_green = 	g <= 1.0 ? (int) Math.floor(g * 255) : 255;
        int c_blue = 	b <= 1.0 ? (int) Math.floor(b * 255) : 255;
        int pixelColor = getColorInteger(c_red,c_green,c_blue,c_alpha);
        return pixelColor;
	}
    //Do NOT modify this function    
    public void resetImage(){
    	for (int j = 0; j < image.getHeight(); j++) {
	        for (int i = 0; i < image.getWidth(); i++) {
	            image.setRGB(i, j, 0);
	        }
	    }
    }
    //Do NOT modify this function
    void getLightVector(double[] lightVector, double[] halfVector, double[] viewVec){
    	VectorMath.setVector(lightVector, -viewVec[0], -viewVec[1], -viewVec[2]);
    	for (int i=0; i<3; i++) {
            halfVector[i] = -viewVec[i] + lightVector[i];
        }
        double l = VectorMath.length(halfVector);
        for (int i=0; i<3; i++) {
            halfVector[i] /= l;
        }
    }
    //Do NOT modify this function
    void computeIncrements(double[] increments, double[] viewVec, double sampleStep) {
    	VectorMath.setVector(increments, viewVec[0] * sampleStep, viewVec[1] * sampleStep, viewVec[2] * sampleStep);
    }
    
    //used by the slicer
    //Do NOT modify this function
    void getPlane(double viewVec[], double uVec[], double vVec[]) {
		VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
	    VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
	    VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
	}
    
    //used by the slicer	
    //Do NOT modify this function
	void computeVolumeCenter(double volumeCenter[]) {
		VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
	}
	
    //used by the slicer
    //Do NOT modify this function
	void computePixelCoordinates(double pixelCoord[], double volumeCenter[], double uVec[], double vVec[], int i, int j) {
		int imageCenter = image.getWidth()/2;
		pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
	}
    //Do NOT modify this function
	void computePixelCoordinatesBehind(double pixelCoord[], double viewVec[], double uVec[], double vVec[], int i, int j) {
		int imageCenter = image.getWidth()/2;
		pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) - viewVec[0] * imageCenter + volume.getDimX() / 2.0;
        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) - viewVec[1] * imageCenter + volume.getDimY() / 2.0;
        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) - viewVec[2] * imageCenter + volume.getDimZ() / 2.0;
	}
	
    //Do NOT modify this function
    public int getColorInteger(int c_red, int c_green, int c_blue, int c_alpha) {
    	int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
    	return pixelColor;
    } 
    //Do NOT modify this function
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }
    //Do NOT modify this function
    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFunc.setTestFunc();
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }
    //Do NOT modify this function
    public RaycastRendererPanel getPanel() {
        return panel;
    }

    //Do NOT modify this function
    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    //Do NOT modify this function
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    //Do NOT modify this function
    public void setShadingMode(boolean mode) {
        shadingMode = mode;
        changed();
    }
    //Do NOT modify this function
    public void setMIPMode() {
        setMode(false, true, false, false);
    }
    //Do NOT modify this function
    public void setSlicerMode() {
        setMode(true, false, false, false);
    }
    //Do NOT modify this function
    public void setCompositingMode() {
        setMode(false, false, true, false);
    }
    //Do NOT modify this function
    public void setTF2DMode() {
        setMode(false, false, false, true);
    }
    //Do NOT modify this function
    private void setMode(boolean slicer, boolean mip, boolean composite, boolean tf2d) {
        slicerMode = slicer;
        mipMode = mip;
        compositingMode = composite;
        tf2dMode = tf2d;        
        changed();
    }
    //Do NOT modify this function
    private boolean intersectLinePlane(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection) {

        double[] tmp = new double[3];

        for (int i = 0; i < 3; i++) {
            tmp[i] = plane_pos[i] - line_pos[i];
        }

        double denom = VectorMath.dotproduct(line_dir, plane_normal);
        if (Math.abs(denom) < 1.0e-8) {
            return false;
        }

        double t = VectorMath.dotproduct(tmp, plane_normal) / denom;

        for (int i = 0; i < 3; i++) {
            intersection[i] = line_pos[i] + t * line_dir[i];
        }

        return true;
    }
    //Do NOT modify this function
    private boolean validIntersection(double[] intersection, double xb, double xe, double yb,
            double ye, double zb, double ze) {

        return (((xb - 0.5) <= intersection[0]) && (intersection[0] <= (xe + 0.5))
                && ((yb - 0.5) <= intersection[1]) && (intersection[1] <= (ye + 0.5))
                && ((zb - 0.5) <= intersection[2]) && (intersection[2] <= (ze + 0.5)));

    }
    //Do NOT modify this function
    private void intersectFace(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection,
            double[] entryPoint, double[] exitPoint) {

        boolean intersect = intersectLinePlane(plane_pos, plane_normal, line_pos, line_dir,
                intersection);
        if (intersect) {

            double xpos0 = 0;
            double xpos1 = volume.getDimX();
            double ypos0 = 0;
            double ypos1 = volume.getDimY();
            double zpos0 = 0;
            double zpos1 = volume.getDimZ();

            if (validIntersection(intersection, xpos0, xpos1, ypos0, ypos1,
                    zpos0, zpos1)) {
                if (VectorMath.dotproduct(line_dir, plane_normal) > 0) {
                    entryPoint[0] = intersection[0];
                    entryPoint[1] = intersection[1];
                    entryPoint[2] = intersection[2];
                } else {
                    exitPoint[0] = intersection[0];
                    exitPoint[1] = intersection[1];
                    exitPoint[2] = intersection[2];
                }
            }
        }
    }
    //Do NOT modify this function
    public double computeLevoyOpacity(double material_value, double material_r,
            double voxelValue, double gradMagnitude) {

        double opacity = 0.0;

        if (gradMagnitude == 0.0 && voxelValue == material_value) {
            opacity = 1.0;
        } else if (gradMagnitude > 0.0 && voxelValue - material_r * gradMagnitude <= material_value
                && material_value <= voxelValue + material_r * gradMagnitude) {
            //opacity = 1.0 - (1.0/material_r) * Math.abs((material_value - voxelValue)/gradMagnitude);
            opacity = 1.0 - Math.abs((material_value - voxelValue) / (gradMagnitude * material_r));
        }

        return opacity;
    }
    //Do NOT modify this function
    public double computeLevoy(int f_l, int f_h, double voxelValue, double gradMagnitude) {
        double opacity = 0.0;
        
        if (f_l <= voxelValue && voxelValue <= f_h) {
            TFColor c_l = tFunc.getColor(f_l);
            TFColor c_h = tFunc.getColor(f_h);
            opacity = gradMagnitude * ((c_h.a * (voxelValue - f_l)/(f_h-f_l)) + (c_l.a * (f_h - voxelValue)/(f_h-f_l)));
        }
        
        return opacity;
    }
    //Do NOT modify this function
    private TFColor computePhongShading(TFColor voxel_color, VoxelGradient gradient, double[] lightVector,
            double[] halfVector) {

        double diffuse_coefficient = 0.7;
        double ambient_coefficient = 0.1;
        double specular_coefficient = 0.2;
        double specular_power = 10;

        double[] grad = new double[3];
        VectorMath.setVector(grad, gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag);

        double diffuse = VectorMath.dotproduct(grad, lightVector);
        
        TFColor color = new TFColor(voxel_color.r, voxel_color.g, voxel_color.b, voxel_color.a);
        
        if (diffuse > 0) {
            color.r = voxel_color.r * diffuse * diffuse_coefficient + ambient_coefficient;
            color.g = voxel_color.g * diffuse * diffuse_coefficient + ambient_coefficient;
            color.b = voxel_color.b * diffuse * diffuse_coefficient + ambient_coefficient;
        }
        double specular = VectorMath.dotproduct(grad, halfVector);
        if (specular > 0) {
            color.r += specular_coefficient * Math.pow(specular, specular_power);
            color.g += specular_coefficient * Math.pow(specular, specular_power);
            color.b += specular_coefficient * Math.pow(specular, specular_power);
        }
        color.r = color.r > 1.0 ? 1.0 : color.r;
        color.g = color.g > 1.0 ? 1.0 : color.g;
        color.b = color.b > 1.0 ? 1.0 : color.b;
        
        return color;
    }
    
    
    //Do NOT modify this function
    void computeEntryAndExit(double[] p, double[] viewVec, double[] entryPoint, double[] exitPoint) {

        for (int i = 0; i < 3; i++) {
            entryPoint[i] = -1;
            exitPoint[i] = -1;
        }

        double[] plane_pos = new double[3];
        double[] plane_normal = new double[3];
        double[] intersection = new double[3];

        VectorMath.setVector(plane_pos, volume.getDimX(), 0, 0);
        VectorMath.setVector(plane_normal, 1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, -1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, volume.getDimY(), 0);
        VectorMath.setVector(plane_normal, 0, 1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, -1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, volume.getDimZ());
        VectorMath.setVector(plane_normal, 0, 0, 1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, 0, -1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

    }

    //Do NOT modify this function
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }
    //Do NOT modify this function
    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        if (slicerMode) {
            slicer(viewMatrix);    
        } else {
            raycast(viewMatrix);
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];
    //Do NOT modify this function
    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
