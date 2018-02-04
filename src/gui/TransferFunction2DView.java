/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import volvis.TransferFunction2D;

/**
 *
 * @author michel
 *  Modified by Anna Vilanova
 */
public class TransferFunction2DView extends javax.swing.JPanel {

    TransferFunction2DEditor ed;
    private final int DOTSIZE = 8;
    public Ellipse2D.Double baseControlPoint, radiusControlPoint;
    boolean selectedBaseControlPoint, selectedRadiusControlPoint;
    
    public Ellipse2D.Double minGradControlPoint, maxGradControlPoint;
    boolean selectedMinGradControlPoint, selectedMaxGradControlPoint;
    int minStart = this.getHeight();
    int maxStart = 0;
    
    /**
     * Creates new form TransferFunction2DView
     * @param ed
     */
    public TransferFunction2DView(TransferFunction2DEditor ed) {
        initComponents();
        
        this.ed = ed;
        selectedBaseControlPoint = false;
        selectedRadiusControlPoint = false;
        selectedMinGradControlPoint = false;
        selectedMaxGradControlPoint = false;
        addMouseMotionListener(new TriangleWidgetHandler());
        addMouseListener(new SelectionHandler());
    }
    
    @Override
    public void paintComponent(Graphics g) {

        Graphics2D g2 = (Graphics2D) g;

        int w = this.getWidth();
        int h = this.getHeight();
        g2.setColor(Color.white);
        g2.fillRect(0, 0, w, h);
        
        double maxHistoMagnitude = ed.histogram[0];
        for (int i = 0; i < ed.histogram.length; i++) {
            maxHistoMagnitude = ed.histogram[i] > maxHistoMagnitude ? ed.histogram[i] : maxHistoMagnitude;
        }
        
        double binWidth = (double) w / (double) ed.xbins;
        double binHeight = (double) h / (double) ed.ybins;
        maxHistoMagnitude = Math.log(maxHistoMagnitude);
        
        
        for (int y = 0; y < ed.ybins; y++) {
            for (int x = 0; x < ed.xbins; x++) {
                if (ed.histogram[y * ed.xbins + x] > 0) {
                    int intensity = (int) Math.floor(255 * (1.0 - Math.log(ed.histogram[y * ed.xbins + x]) / maxHistoMagnitude));
                    g2.setColor(new Color(intensity, intensity, intensity));
                    g2.fill(new Rectangle2D.Double(x * binWidth, h - (y * binHeight), binWidth, binHeight));
                }
            }
        }
        
        int ypos = h;
        int xpos = (int) (ed.tf2D.baseIntensity * binWidth);
        g2.setColor(Color.black);
        baseControlPoint = new Ellipse2D.Double(xpos - DOTSIZE / 2, ypos - DOTSIZE, DOTSIZE, DOTSIZE);
        g2.fill(baseControlPoint);

        g2.drawLine(xpos, ypos, xpos - (int) (ed.tf2D.radius * binWidth * ed.maxGradientMagnitude), 0);
        g2.drawLine(xpos, ypos, xpos + (int) (ed.tf2D.radius * binWidth * ed.maxGradientMagnitude), 0);
        radiusControlPoint = new Ellipse2D.Double(xpos + (ed.tf2D.radius * binWidth * ed.maxGradientMagnitude) - DOTSIZE / 2,  0, DOTSIZE, DOTSIZE);
        g2.fill(radiusControlPoint);
        
        if(minStart == 0){
            minStart = h;
        }
        minGradControlPoint = new Ellipse2D.Double(0, minStart - DOTSIZE / 2, DOTSIZE, DOTSIZE);
        maxGradControlPoint = new Ellipse2D.Double(0, maxStart - DOTSIZE / 2 , DOTSIZE, DOTSIZE);
        g2.drawLine(0,minStart,w,minStart);
        g2.drawLine(0,maxStart,w,maxStart);
        g2.fill(minGradControlPoint);
        g2.fill(maxGradControlPoint);
    }
    
    
    private class TriangleWidgetHandler extends MouseMotionAdapter {

        @Override
        public void mouseMoved(MouseEvent e) {
            if (baseControlPoint.contains(e.getPoint()) || 
                    radiusControlPoint.contains(e.getPoint()) || 
                    minGradControlPoint.contains(e.getPoint()) ||
                    maxGradControlPoint.contains(e.getPoint())) {
                setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            } else {
                setCursor(Cursor.getDefaultCursor());
            }
        }
        
        @Override
        public void mouseDragged(MouseEvent e) {
            if (selectedBaseControlPoint || selectedRadiusControlPoint ||
                    selectedMinGradControlPoint ||
                    selectedMaxGradControlPoint) {
                Point dragEnd = e.getPoint();
                
                if (selectedBaseControlPoint) {
                    // restrain to horizontal movement
                    dragEnd.setLocation(dragEnd.x, baseControlPoint.getCenterY());
                } else if (selectedRadiusControlPoint) {
                    // restrain to horizontal movement and avoid radius getting 0
                    dragEnd.setLocation(dragEnd.x, radiusControlPoint.getCenterY());
                    if (dragEnd.x - baseControlPoint.getCenterX() <= 0) {
                        dragEnd.x = (int) (baseControlPoint.getCenterX() + 1);
                    }
                }else if (selectedMinGradControlPoint){
                    dragEnd.setLocation(0, dragEnd.y);
                }else if (selectedMaxGradControlPoint){
                    dragEnd.setLocation(0, dragEnd.y);
                }
                if (dragEnd.x < 0) {
                    dragEnd.x = 0;
                }
                if (dragEnd.x >= getWidth()) {
                    dragEnd.x = getWidth() - 1;
                }
                double w = getWidth();
                double h = getHeight();
                double binWidth = (double) w / (double) ed.xbins;
                if (selectedBaseControlPoint) {
                    ed.tf2D.baseIntensity = (short) (dragEnd.x / binWidth);
                } else if (selectedRadiusControlPoint) {
                    ed.tf2D.radius = (dragEnd.x - (ed.tf2D.baseIntensity * binWidth))/(binWidth*ed.maxGradientMagnitude);
                } else if (selectedMaxGradControlPoint){
                    if(dragEnd.y < minStart){
                        maxStart = dragEnd.y;
                    }else{
                        maxStart = minStart - 1;
                    }
                    ed.tf2D.max = (1 - maxStart/h)*ed.maxGradientMagnitude; 
                } else if (selectedMinGradControlPoint){
                    if(dragEnd.y > maxStart){
                        minStart = dragEnd.y;
                    }else{
                        minStart = maxStart + 1;
                    }
                    ed.tf2D.min = (1 - minStart/h)*ed.maxGradientMagnitude; 
                }
                ed.setSelectedInfo();
                
                repaint();
            } 
        }

    }
    
    
    private class SelectionHandler extends MouseAdapter {
        @Override
        public void mousePressed(MouseEvent e) {
            if (baseControlPoint.contains(e.getPoint())) {
                selectedBaseControlPoint = true;
            } else if (radiusControlPoint.contains(e.getPoint())) {
                selectedRadiusControlPoint = true;
            } else if (minGradControlPoint.contains(e.getPoint())) {
                selectedMinGradControlPoint = true;
            } else if (maxGradControlPoint.contains(e.getPoint())) {
                selectedMaxGradControlPoint = true;
            } else {
                selectedRadiusControlPoint = false;
                selectedBaseControlPoint = false;
                selectedMaxGradControlPoint = false;
                selectedMinGradControlPoint = false;
            }
        }
        
        @Override
        public void mouseReleased(MouseEvent e) {
            selectedRadiusControlPoint = false;
            selectedBaseControlPoint = false;
            selectedMaxGradControlPoint = false;
            selectedMinGradControlPoint = false;
            ed.changed();
            repaint();
        }
    }
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents


    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}
