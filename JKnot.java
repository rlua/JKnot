/**
 * @author Rhonald Lua
 * @version 1.0 started 01/26/2002
 * last modified 06/20/2002
 */

// Assumptions: polymer fills all lattice, dimension side^3 (recently relaxed),
// each polymer joint only connected to two nearest neighbors, one segment each.

import java.awt.*;
import java.applet.*;
import java.awt.event.*;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.Enumeration;
import java.math.BigInteger;

// Node (or vertex) always connected to two other nodes, one segment each
class JNode
{
	public int prevsegindex,	// segments of path sequentially indexed
		nextsegindex;
	public char	prevsegdir,	// direction to previous node	'F','B','U','D','L','R'
		nextsegdir;	// direction to next node

	JNode()
	{
		prevsegindex=-1;
		nextsegindex=-1;
		prevsegdir='N';
		nextsegdir='N';
	}
}

class JCross
{
	public boolean type;	// over (true) or under (false)
	public int over,under;	// segments that make up the crossing
	public double x,y;	// coordinates of intersection in plane projection
	public boolean xytype;	// true if the vertical segment is over the horizontal segment
	public char underdir;	// Used to determine underpass type. see A.V. Vologodskii, et al, Sov. Phys.-JETP, Vol. 39, 1974, p. 1059
	public char overdir;
	public int partner;
	public int underpassnum;
	public int generatornum;
	public byte crosstype;	// Used to calculate writhing number. see Deguchi and Tsurusaki, Phys. Letters A 174, 1993, 29-37

	public Object partnerObj;	// should be a reference to the JCross object that serves as a partner, i.e. that intersects at the same point

	JCross(boolean b, int o, int u, double xo, double yo, boolean b2, char odir, char udir)
	{
		type=b;
		over=o;
		under=u;

		x=xo;
		y=yo;
		xytype=b2;

		underdir=udir;
		overdir=odir;

		partner=-1;
		underpassnum=0;
		generatornum=0;

		crosstype=1;

		partnerObj=null;
	}
}

class JPolymer
{
	public int[] x,y,z;
	public int size;	//number of points; (lattice size)^3
	public int side;

	JPolymer()
	{
		size=0;
	}

	public boolean Parse(String s, int c)
	{
		StringTokenizer st=new StringTokenizer(s,"\r\n");
		size=st.countTokens();
		if(size>0)
		{
			x=new int[size];
			y=new int[size];
			z=new int[size];
			//side=(int)Math.ceil(Math.pow(size,1.0/3.0));
		}
		int i=0,xyzmax=0;
		while(st.hasMoreTokens())
		{
			StringTokenizer st2=new StringTokenizer(st.nextToken()," ");
			if(st2.countTokens()<3)
				return false;
			try
			{
				// try different projections (xy,zx,yz)
				if(c==0)
				{
					x[i]=Integer.parseInt(st2.nextToken());
					y[i]=Integer.parseInt(st2.nextToken());
					z[i]=Integer.parseInt(st2.nextToken());
				}
				else if(c==1)
				{
					y[i]=Integer.parseInt(st2.nextToken());
					z[i]=Integer.parseInt(st2.nextToken());
					x[i]=Integer.parseInt(st2.nextToken());
				}
				else if(c==2)
				{
					z[i]=Integer.parseInt(st2.nextToken());
					x[i]=Integer.parseInt(st2.nextToken());
					y[i]=Integer.parseInt(st2.nextToken());
				}
				if(xyzmax<x[i])
					xyzmax=x[i];
				if(xyzmax<y[i])
					xyzmax=y[i];
				if(xyzmax<z[i])
					xyzmax=z[i];
				i++;
			}
			catch(Exception ex)
			{
				return false;
			}
		}

		side=xyzmax+1;	//allow for non-cubic conformation

		return true;
	}
}

class Coords
{
	public double x,y,z;

	Coords()
	{

	}

	Coords(double xo, double yo, double zo)
	{
		x=xo;
		y=yo;
		z=zo;
	}
}

// A 3D perspective of polymer (left panel)
class J3D extends Canvas implements MouseMotionListener
{
	public JPolymer polyref;

	double imagesize;
	double objectsize;

	boolean WithEyetoScreen;
	boolean ShowAxes;

	//phi and theta convention interchanged?
	final double PI = 3.14159;
	final int TABLESIZE = 60;

	// angle indices
	double rho;
	int phi;
	int theta;

	double[] sintab=new double[TABLESIZE];
	double[] costab=new double[TABLESIZE];

	J3D(int length)
	{
		//fix view point
		objectsize=length/2;	// characteristic length
		imagesize=3*objectsize/4;

		WithEyetoScreen=true;
		ShowAxes=false;

		rho=10*objectsize;
		phi=0;
		theta=0;

		inittabs();

		addMouseMotionListener(this);
	}

	//initialize sine and cosine tables, for fast computation
	void inittabs()
	{
		for(int i=0;i<TABLESIZE;i++)
		{
			sintab[i]=Math.sin(2*PI*i/TABLESIZE);	//increments of 6 degrees
			costab[i]=Math.cos(2*PI*i/TABLESIZE);
		}
	}

	// viewing transformation, from world coords to eye coords
	Coords WorldtoEye( Coords w, double rho, int phi, int theta )
	{
		Coords e=new Coords();
		e.x = -w.x*sintab[phi]+w.y*costab[phi];
		e.y = -w.x*costab[theta]*costab[phi]-w.y*costab[theta]*sintab[phi]+w.z*sintab[theta];
		e.z = -w.x*sintab[theta]*costab[phi]-w.y*sintab[theta]*sintab[phi]-w.z*costab[theta]+rho;
		return e;
	}

	// perspective transformation, eye to screen coords
	Coords EyetoScreen( Coords e , double rho, double imsize, double objsize )
	{
		double d=rho*imsize/objsize;
		Coords s=new Coords();
		s.x = Math.ceil(d*e.x/e.z);
		s.y = Math.ceil(d*e.y/e.z);
		return s;
	}

	void ThreeDLine(Coords p1, Coords p2, Graphics g)
	{
		Coords Eye=new Coords();
		Coords Screen=new Coords();
		int x1,y1,x2,y2;

		//World=p1
		Eye = WorldtoEye( p1, rho, phi, theta );
		if(WithEyetoScreen)
		{
			Screen = EyetoScreen( Eye, rho, imagesize, objectsize );
			x1=(int)Screen.x;
			y1=(int)Screen.y;
		}
		else
		{
			x1=(int)Eye.x;
			y1=(int)Eye.y;
		}

		//World=p2
		Eye = WorldtoEye( p2, rho, phi, theta );
		if(WithEyetoScreen)
		{
			Screen = EyetoScreen( Eye, rho, imagesize, objectsize );
			x2=(int)Screen.x;
			y2=(int)Screen.y;
		}
		else
		{
			x2=(int)Eye.x;
			y2=(int)Eye.y;
		}

		// shift origin here
		g.drawLine(x1+this.getSize().width/2,y1+this.getSize().height/2,x2+this.getSize().width/2,y2+this.getSize().height/2);
	}

	public void HandleEvent(String label)
	{
		if(label.equals("up"))
		{
			theta+=1;
			theta%=TABLESIZE;
			repaint();
		}
		else if(label.equals("down"))
		{
			theta-=1;
			theta+=TABLESIZE;
			theta%=TABLESIZE;
			repaint();
		}
		else if(label.equals("left"))
		{
			phi-=1;
			phi+=TABLESIZE;
			phi%=TABLESIZE;
			repaint();
		}
		else if(label.equals("right"))
		{
			phi+=1;
			phi%=TABLESIZE;
			repaint();
		}
		else if(label.equals("zoom in"))
		{
			rho-=objectsize;
			repaint();
		}
		else if(label.equals("zoom out"))
		{
			rho+=objectsize;
			repaint();
		}
		else if(label.equals("isom"))
		{
			WithEyetoScreen=!WithEyetoScreen;
			repaint();
		}
		else if(label.equals("axes"))
		{
			ShowAxes=!ShowAxes;
			repaint();
		}

	}

	int Xprev=0,Yprev=0;
	public void mouseDragged(MouseEvent evt)
	{
		// called when the mouse is moved while
		// a button is pressed down
		int x=evt.getX();
		int y=evt.getY();

		if(y-Yprev>1)
		{
			HandleEvent("up");
		}
		else if(y-Yprev<-1)
		{
			HandleEvent("down");
		}
		Xprev=x;
		Yprev=y;
	}

	public void mouseMoved(MouseEvent evt)
	{
		// called when the mouse is moved while
		// no buttons are pressed down
		Xprev=evt.getX();
		Yprev=evt.getY();
	}

	public void Reset()
	{
		phi=0;
		theta=0;
	}

	public void paint(Graphics g)
	{
		if(ShowAxes)
		{
			//Color oldColor;
			//oldColor=g.getColor();
			g.setColor(Color.gray);
			ThreeDLine(new Coords(objectsize,0,0),new Coords(-objectsize,0,0),g);
			g.setColor(Color.green);
			ThreeDLine(new Coords(0,objectsize,0),new Coords(0,-objectsize,0),g);
			g.setColor(Color.blue);
			ThreeDLine(new Coords(0,0,objectsize),new Coords(0,0,-objectsize),g);
			//g.setColor(oldColor);
		}

		Coords c1=new Coords(0,0,0);
		Coords c2=new Coords(0,0,0);
		int i=0;
		if(polyref.size>1)
		{
			g.setColor(Color.yellow);	// different color for starting segment
			c1.x=(objectsize*polyref.x[i])/polyref.side;
			c1.y=(objectsize*polyref.y[i])/polyref.side;
			c1.z=(objectsize*polyref.z[i])/polyref.side;
			c2.x=(objectsize*polyref.x[i+1])/polyref.side;
			c2.y=(objectsize*polyref.y[i+1])/polyref.side;
			c2.z=(objectsize*polyref.z[i+1])/polyref.side;

			ThreeDLine(c1,c2,g);
			i++;
		}
		g.setColor(Color.red);
		for(i=1;i<polyref.size-1;i++)
		{
			c1.x=(objectsize*polyref.x[i])/polyref.side;
			c1.y=(objectsize*polyref.y[i])/polyref.side;
			c1.z=(objectsize*polyref.z[i])/polyref.side;
			c2.x=(objectsize*polyref.x[i+1])/polyref.side;
			c2.y=(objectsize*polyref.y[i+1])/polyref.side;
			c2.z=(objectsize*polyref.z[i+1])/polyref.side;

			ThreeDLine(c1,c2,g);
		}
	}


}

// The plane projection (right panel)
class J2D extends Canvas
{
	public JPolymer polyref;
	public Vector crossref;

	public void paint(Graphics g)
	{
		int w=this.getSize().width-10;
		int h=this.getSize().height-10;
		int i=0,x1,y1,x2,y2;
		if(polyref.size>1)
		{
			g.setColor(Color.red);	// different color for starting segment
			x1=4+(w*polyref.x[i]+(w*polyref.z[i])/polyref.side)/polyref.side;
			y1=4+(h*polyref.y[i]+(h*polyref.z[i])/polyref.side)/polyref.side;
			x2=4+(w*polyref.x[i+1]+(w*polyref.z[i+1])/polyref.side)/polyref.side;
			y2=4+(h*polyref.y[i+1]+(h*polyref.z[i+1])/polyref.side)/polyref.side;
			g.drawLine(x1,y1,x2,y2);
			i++;
		}
		g.setColor(Color.yellow);
		// Use transformation suggested by Alex: x -> x + z/n, y -> y + z/n
		for(i=1;i<polyref.size-1;i++)
		{
			x1=4+(w*polyref.x[i]+(w*polyref.z[i])/polyref.side)/polyref.side;	// already calculated before? x,y
			y1=4+(h*polyref.y[i]+(h*polyref.z[i])/polyref.side)/polyref.side;
			x2=4+(w*polyref.x[i+1]+(w*polyref.z[i+1])/polyref.side)/polyref.side;
			y2=4+(h*polyref.y[i+1]+(h*polyref.z[i+1])/polyref.side)/polyref.side;
			g.drawLine(x1,y1,x2,y2);
		}

		if(polyref.size>1)
		{
			for(i=0;i<crossref.size();i++)	// double drawing?
			{
				JCross jtmp=(JCross)crossref.elementAt(i);
				if(jtmp.type)	continue;
				if(jtmp.xytype)
				{
					g.setColor(Color.white);	// vertical segment over horizontal segment
				}
				else
				{
					g.setColor(Color.black);	// vertical segment under horizontal segment
				}
				g.drawArc(2+(int)(w*jtmp.x/polyref.side),2+(int)(h*jtmp.y/polyref.side),4,4,0,360);
			}
		}
		else if(polyref.size==-1)	// draw chord diagram, assumes partners filled up by Dowker
		{
			int cx=this.getSize().width;
			int cy=this.getSize().height;
			int cc=(cx<cy?cx:cy);
			int R=cc-10;
			g.setColor(Color.yellow);
			g.drawArc(5,5,R,R,0,360);
			R/=2;
			for(i=0;i<crossref.size();i++)
			{
				JCross jtmp=(JCross)crossref.elementAt(i);
				if(jtmp.type)	continue;
				x1=5+R+(int)(R*Math.cos(i*2*3.14159/crossref.size()));
				y1=5+R+(int)(R*Math.sin(i*2*3.14159/crossref.size()));
				g.setColor(Color.red);
				g.drawArc(x1-4,y1-4,8,8,0,360);	// draw circle around underpass crossing
				if(jtmp.crosstype>0)	// writhe; (+), (-)
				{
					g.setColor(Color.orange);
				}
				else
				{
					g.setColor(Color.black);
				}
				x2=5+R+(int)(R*Math.cos(jtmp.partner*2*3.14159/crossref.size()));
				y2=5+R+(int)(R*Math.sin(jtmp.partner*2*3.14159/crossref.size()));
				g.drawLine(x1,y1,x2,y2);
			}
		}
	}
}

/**
 * The applet housing
 */
public class JKnot extends Applet
{
	J3D cv3D;
	J2D cv2D=new J2D();
	Button btInput=new Button("Process Input");
	Button btUp=new Button("up");
	Button btDown=new Button("down");
	Button btLeft=new Button("left");
	Button btRight=new Button("right");
	Button btClear=new Button("Clear Output");
	TextArea taInput=new TextArea(),taOutput=new TextArea();
	TextField tfAlexVar=new TextField("-1");
	Checkbox cbDet=new Checkbox("Determinant, t =",true);
	Checkbox cbReduce=new Checkbox("Reduce",true);
	Choice chInput;
	Button btClearInput=new Button("Clear Input");
	JPolymer Polymer=new JPolymer();

	JNode[][][] Nodes;	// Nodes[x][y][z]
	Vector Crossing=new Vector();

	// Fill up Nodes and Crossing;  Used to generate Dowker notation.  Assume knot fills all of lattice
	// Path: ... size-1 -> 0 -> 1 -> 2 ...
	void ParsePolymer() throws Exception
	{
		int i,x=0,y=0,z=0;
		//Hashtable[] Cross;	// map of path segments to crossing type; 'U' (1 under 2), 'O' (1 over 2)
		//Cross=new Hashtable[Polymer.size];	// phase out?
		Nodes=new JNode[Polymer.side][Polymer.side][Polymer.side];
		for(x=0;x<Polymer.side;x++)	// initialize to allow for non-cubic conformations
		{
			for(y=0;y<Polymer.side;y++)
			{
				for(z=0;z<Polymer.side;z++)
				{
					Nodes[x][y][z]=new JNode();
				}
			}
		}
		Crossing.removeAllElements();	// clear Crossing first
		for(i=0;i<Polymer.size;i++)
		{
			//Cross[i]=new Hashtable();
			int tmp;
			char dtmp;
			x=Polymer.x[i];	// current node
			y=Polymer.y[i];
			z=Polymer.z[i];

			if(i==0)
			{
				tmp=Polymer.size-1;
			}
			else
			{
				tmp=i-1;
			}
			//Nodes[x][y][z]=new JNode();	// ah, yes...
			Nodes[x][y][z].prevsegindex=tmp;
			dtmp=SegDir(x,y,z,Polymer.x[tmp],Polymer.y[tmp],Polymer.z[tmp]);
			if(dtmp=='N')
			{
				taOutput.append("Error parsing polymer:\n Incorrect segment encountered\n");
				return;
			}
			Nodes[x][y][z].prevsegdir=dtmp;
			if(i==Polymer.size-1)
			{
				tmp=0;
			}
			else
			{
				tmp=i+1;
			}
			Nodes[x][y][z].nextsegindex=i;
			dtmp=SegDir(x,y,z,Polymer.x[tmp],Polymer.y[tmp],Polymer.z[tmp]);
			if(dtmp=='N')
			{
				taOutput.append("Error parsing polymer:\n Incorrect segment encountered\n");
				return;
			}
			Nodes[x][y][z].nextsegdir=dtmp;
		}

		int ztmp,tmpindex;
		double pside=Polymer.side;
		for(i=0;i<Polymer.size;i++)
		{
			// try different projections (xy,zx,yz)
			x=Polymer.x[i];	// current node
			y=Polymer.y[i];
			z=Polymer.z[i];

			// Only horizontals and verticals get crossed in the projection
			if(Nodes[x][y][z].nextsegdir=='F')
			{
				// checking sequence important?
				// for the smaller x endpoint of the segment check for y direction of upper z
				for(ztmp=z+1;ztmp<Polymer.side;ztmp++)
				{
					if(Nodes[x][y][ztmp].nextsegdir=='L')
					{
						tmpindex=Nodes[x][y][ztmp].nextsegindex;
						//Cross[i].put("s"+tmpindex,"U");
						//Cross[tmpindex].put("s"+i,"O");
/*1*/						Crossing.addElement(new JCross(false,tmpindex,i,x+ztmp/pside,y+z/pside,true,'L','F'));	// (x-vertical,y-horizontal)
					}
					if(Nodes[x][y][ztmp].prevsegdir=='L')
					{
						tmpindex=Nodes[x][y][ztmp].prevsegindex;
						//Cross[i].put("s"+tmpindex,"U");
						//Cross[tmpindex].put("s"+i,"O");
/*2*/						Crossing.addElement(new JCross(false,tmpindex,i,x+ztmp/pside,y+z/pside,true,'R','F'));
					}
				}
				// for the larger x endpoint of the segment check for y direction of lower z
				for(ztmp=0;ztmp<z;ztmp++)
				{
					if(Nodes[x+1][y][ztmp].nextsegdir=='R')
					{
						tmpindex=Nodes[x+1][y][ztmp].nextsegindex;
						//Cross[i].put("s"+tmpindex,"O");
						//Cross[tmpindex].put("s"+i,"U");
/*3*/						Crossing.addElement(new JCross(true,i,tmpindex,x+1+ztmp/pside,y+z/pside,false,'F','R'));
					}
					if(Nodes[x+1][y][ztmp].prevsegdir=='R')
					{
						tmpindex=Nodes[x+1][y][ztmp].prevsegindex;
						//Cross[i].put("s"+tmpindex,"O");
						//Cross[tmpindex].put("s"+i,"U");
/*4*/						Crossing.addElement(new JCross(true,i,tmpindex,x+1+ztmp/pside,y+z/pside,false,'F','L'));
					}
				}
			}
			if(Nodes[x][y][z].nextsegdir=='B')
			{
				for(ztmp=z-1;ztmp>=0;ztmp--)
				{
					if(Nodes[x][y][ztmp].nextsegdir=='R')
					{
						tmpindex=Nodes[x][y][ztmp].nextsegindex;
						//Cross[i].put("s"+tmpindex,"O");
						//Cross[tmpindex].put("s"+i,"U");
/*5*/						Crossing.addElement(new JCross(true,i,tmpindex,x+ztmp/pside,y+z/pside,false,'B','R'));
					}
					if(Nodes[x][y][ztmp].prevsegdir=='R')
					{
						tmpindex=Nodes[x][y][ztmp].prevsegindex;
						//Cross[i].put("s"+tmpindex,"O");
						//Cross[tmpindex].put("s"+i,"U");
/*6*/						Crossing.addElement(new JCross(true,i,tmpindex,x+ztmp/pside,y+z/pside,false,'B','L'));
					}
				}
				for(ztmp=Polymer.side-1;ztmp>z;ztmp--)
				{
					if(Nodes[x-1][y][ztmp].nextsegdir=='L')
					{
						tmpindex=Nodes[x-1][y][ztmp].nextsegindex;
						//Cross[i].put("s"+tmpindex,"U");
						//Cross[tmpindex].put("s"+i,"O");
/*7*/						Crossing.addElement(new JCross(false,tmpindex,i,x-1+ztmp/pside,y+z/pside,true,'L','B'));
					}
					if(Nodes[x-1][y][ztmp].prevsegdir=='L')
					{
						tmpindex=Nodes[x-1][y][ztmp].prevsegindex;
						//Cross[i].put("s"+tmpindex,"U");
						//Cross[tmpindex].put("s"+i,"O");
/*8*/						Crossing.addElement(new JCross(false,tmpindex,i,x-1+ztmp/pside,y+z/pside,true,'R','B'));
					}
				}
			}
			if(Nodes[x][y][z].nextsegdir=='R')
			{
				// for the smaller y endpoint of the segment check for x direction of upper z
				for(ztmp=z+1;ztmp<Polymer.side;ztmp++)
				{
					if(Nodes[x][y][ztmp].nextsegdir=='B')
					{
						tmpindex=Nodes[x][y][ztmp].nextsegindex;
						//Cross[i].put("s"+tmpindex,"U");
						//Cross[tmpindex].put("s"+i,"O");
/*9*/						Crossing.addElement(new JCross(false,tmpindex,i,x+z/pside,y+ztmp/pside,false,'B','R'));
					}
					if(Nodes[x][y][ztmp].prevsegdir=='B')
					{
						tmpindex=Nodes[x][y][ztmp].prevsegindex;
						//Cross[i].put("s"+tmpindex,"U");
						//Cross[tmpindex].put("s"+i,"O");
/*10*/						Crossing.addElement(new JCross(false,tmpindex,i,x+z/pside,y+ztmp/pside,false,'F','R'));
					}
				}
				// for the larger y endpoint of the segment check for x direction of lower z
				for(ztmp=0;ztmp<z;ztmp++)
				{
					if(Nodes[x][y+1][ztmp].nextsegdir=='F')
					{
						tmpindex=Nodes[x][y+1][ztmp].nextsegindex;
						//Cross[i].put("s"+tmpindex,"O");
						//Cross[tmpindex].put("s"+i,"U");
/*11*/						Crossing.addElement(new JCross(true,i,tmpindex,x+z/pside,y+1+ztmp/pside,true,'R','F'));
					}
					if(Nodes[x][y+1][ztmp].prevsegdir=='F')
					{
						tmpindex=Nodes[x][y+1][ztmp].prevsegindex;
						//Cross[i].put("s"+tmpindex,"O");
						//Cross[tmpindex].put("s"+i,"U");
/*12*/						Crossing.addElement(new JCross(true,i,tmpindex,x+z/pside,y+1+ztmp/pside,true,'R','B'));
					}
				}
			}
			if(Nodes[x][y][z].nextsegdir=='L')
			{
				// for the larger y endpoint of the segment check for x direction of lower z
				for(ztmp=z-1;ztmp>=0;ztmp--)
				{
					if(Nodes[x][y][ztmp].nextsegdir=='F')
					{
						tmpindex=Nodes[x][y][ztmp].nextsegindex;
						//Cross[i].put("s"+tmpindex,"O");
						//Cross[tmpindex].put("s"+i,"U");
/*13*/						Crossing.addElement(new JCross(true,i,tmpindex,x+z/pside,y+ztmp/pside,true,'L','F'));
					}
					if(Nodes[x][y][ztmp].prevsegdir=='F')
					{
						tmpindex=Nodes[x][y][ztmp].prevsegindex;
						//Cross[i].put("s"+tmpindex,"O");
						//Cross[tmpindex].put("s"+i,"U");
/*14*/						Crossing.addElement(new JCross(true,i,tmpindex,x+z/pside,y+ztmp/pside,true,'L','B'));
					}
				}
				// for the smaller y endpoint of the segment check for x direction of upper z
				for(ztmp=Polymer.side-1;ztmp>z;ztmp--)
				{
					if(Nodes[x][y-1][ztmp].nextsegdir=='B')
					{
						tmpindex=Nodes[x][y-1][ztmp].nextsegindex;
						//Cross[i].put("s"+tmpindex,"U");
						//Cross[tmpindex].put("s"+i,"O");
/*15*/						Crossing.addElement(new JCross(false,tmpindex,i,x+z/pside,y-1+ztmp/pside,false,'B','L'));
					}
					if(Nodes[x][y-1][ztmp].prevsegdir=='B')
					{
						tmpindex=Nodes[x][y-1][ztmp].prevsegindex;
						//Cross[i].put("s"+tmpindex,"U");
						//Cross[tmpindex].put("s"+i,"O");
/*16*/						Crossing.addElement(new JCross(false,tmpindex,i,x+z/pside,y-1+ztmp/pside,false,'F','L'));
					}
				}
			}
		}

		if(Crossing.size()%2!=0)
		{
			taOutput.append("Error parsing polymer:\n Total number of underpasses and overpasses should be even\n");
			return;
		}

		taOutput.append("Number of crossings in projection: "+Crossing.size()/2+"\n");

		taOutput.append(MatchPartnerObjects());

		if(cbReduce.getState())
		{
			ReduceCrossings();
			taOutput.append("Number of crossings after reduction: "+Crossing.size()/2+"\n");
		}

		if(Crossing.size()>1)
		{
			taOutput.append(GenerateDowker());
			taOutput.append(GenerateAlexander());
		}

	}	// ParsePolymer

	String MatchPartnerObjects()
	{
		int i,k;
		for(i=0;i<Crossing.size();i++)
		{
			JCross jtmpi=(JCross)Crossing.elementAt(i);
			// Search for partner of crossing i (forms a pair, one odd, one even)
			for(k=0;k<Crossing.size();k++)
			{
				JCross jtmpk=(JCross)Crossing.elementAt(k);
				if(k!=i)
				{
					if(jtmpk.over==jtmpi.over && jtmpk.under==jtmpi.under)
					{
						//jtmpi.partner=k;
						jtmpi.partnerObj=(Object)jtmpk;
						// assert that jtmpi.type and jtmpk.type should be opposite
						if(jtmpk.type==jtmpi.type)
							return "Error:\n Type should be opposite\n";
						break;
					}
				}
			}
			if(k==Crossing.size())
				return "Error:\n Incomplete partner found.\n";
		}
		return "Partner objects matched\n";
	}

	// Eliminate trivial intersections via Reidemeister moves
	void ReduceCrossings() throws Exception
	{
		boolean found=true,macrofound=true;;
		int i=0;
		while(macrofound)
		{
		macrofound=false;
		found=true;
		while(found)
		{
			found=false;
			// Ignore trivial intersections/loops of the Reidemeister I type.  Except wraparound at endpoints?
			i=0;
			while(Crossing.size()>i+1)
			{
				JCross jtmpi1=(JCross)Crossing.elementAt(i);
				JCross jtmpi2=(JCross)Crossing.elementAt(i+1);
				if(jtmpi1.over==jtmpi2.over && jtmpi1.under==jtmpi2.under)	// look for consecutive crossings
				{
					Crossing.removeElementAt(i);
					Crossing.removeElementAt(i);
					i--;
					if(i<0)
						i++;
				}
				else
				{
					i++;
				}
			}
			// Ignore trivial intersections/loops of the Reidemeister II type.
			i=0;
			while(Crossing.size()>i+1)
			{
				JCross jtmpi1=(JCross)Crossing.elementAt(i);
				JCross jtmpi2=(JCross)Crossing.elementAt(i+1);
				if(jtmpi1.type==jtmpi2.type)
				{
					int j=i+2;
					while(Crossing.size()>j+1)	// look for partner crossings
					{
						JCross jtmpi3=(JCross)Crossing.elementAt(j);
						JCross jtmpi4=(JCross)Crossing.elementAt(j+1);
						if( (jtmpi1.over==jtmpi3.over && jtmpi1.under==jtmpi3.under &&
							jtmpi2.over==jtmpi4.over && jtmpi2.under==jtmpi4.under) ||
							(jtmpi1.over==jtmpi4.over && jtmpi1.under==jtmpi4.under &&
							jtmpi2.over==jtmpi3.over && jtmpi2.under==jtmpi3.under) )
						{
							Crossing.removeElementAt(j);
							Crossing.removeElementAt(j);
							Crossing.removeElementAt(i);
							Crossing.removeElementAt(i);
							i--;
							if(i<0)
								i++;
							found=true;
							break;
						}
						j++;
					}
					if(Crossing.size()<=j+1)
						i++;
				}
				else
				{
					i++;
				}
			}
			// Reduction by Reidemeister III (and Reidemeister I), e.g. a one twist link
			i=0;
			while(Crossing.size()>i+3)
			{
				JCross jtmpi1=(JCross)Crossing.elementAt(i);
				JCross jtmpi4=(JCross)Crossing.elementAt(i+3);
				if(jtmpi1.over==jtmpi4.over && jtmpi1.under==jtmpi4.under)
				{
					JCross jtmpi2=(JCross)Crossing.elementAt(i+1);
					JCross jtmpi3=(JCross)Crossing.elementAt(i+2);
					if(jtmpi1.type==jtmpi2.type && jtmpi1.type!=jtmpi3.type)
					{
						int j=0;
						while(Crossing.size()>j+1)	// also look for and swap the partner crossings
						{
							if(j+1<i || j>i+3)
							{
								JCross jtmpi2p=(JCross)Crossing.elementAt(j);
								JCross jtmpi3p=(JCross)Crossing.elementAt(j+1);
								if( (jtmpi2.over==jtmpi2p.over && jtmpi2.under==jtmpi2p.under &&
									jtmpi3.over==jtmpi3p.over && jtmpi3.under==jtmpi3p.under) ||
									(jtmpi2.over==jtmpi3p.over && jtmpi2.under==jtmpi3p.under &&
									jtmpi3.over==jtmpi2p.over && jtmpi3.under==jtmpi2p.under) )
								{
									JCross jswap23p=jtmpi2p;
									Crossing.setElementAt(jtmpi3p,j);
									Crossing.setElementAt(jswap23p,j+1);

									Crossing.removeElementAt(i+3);
									Crossing.removeElementAt(i);
									found=true;

									break;
								}
							}
							j++;
						}
						if(Crossing.size()<=j+1) i++;
					}
					else
					{
						i++;
					}
				}
				else
				{
					i++;
				}
			}

		// test 'macro' Reidemeister move
		Vector CondemnedCrossings=new Vector();
		int j,k;
		i=0;
		while(Crossing.size()>i+5)
		{
			JCross jtmpi1=(JCross)Crossing.elementAt(i);
			JCross jtmpi2=(JCross)Crossing.elementAt(i+1);
			JCross jtmpi3=null;
			for(j=i+2;j<Crossing.size();j++)
			{
				jtmpi3=(JCross)Crossing.elementAt(j);
				if(jtmpi3.type!=jtmpi2.type)	// 'streak' broken
					break;
				if(jtmpi3.over==jtmpi1.over && jtmpi3.under==jtmpi1.under)	// back to start, 'fish' completed
					break;
			}
			if(j<Crossing.size())
			{
				if(j-i<5)
				{
					i=j;
				}
				else
				{
					if(jtmpi3.over==jtmpi1.over && jtmpi3.under==jtmpi1.under)
					{
						for(k=j-1;k>i;k--)	// remove crossings
						{
							jtmpi3=(JCross)Crossing.elementAt(k);
							CondemnedCrossings.addElement(jtmpi3.partnerObj);
							Crossing.removeElementAt(k);
						}
						for(k=0;k<CondemnedCrossings.size();k++)
						{
							//Crossing.removeElementAt(Crossing.indexOf(CondemnedCrossings.elementAt(k)));
							Crossing.removeElement(CondemnedCrossings.elementAt(k));
						}
						CondemnedCrossings.removeAllElements();
						macrofound=true;
					}
					else
					{
						i++;
					}
				}
			}
			else
			{
				i++;
			}
		}

		} //while(found)
		} //while(macrofound)

	}

	// Idea based on The Knot Book, by Adams
	String GenerateDowker()
	{
		int istart,i,k;
		//int[] partner=new int[Crossing.size()];	// partner[crossing1]=crossing2
		for(istart=0;istart<Crossing.size();istart++)
		{
			JCross jtmp=(JCross)Crossing.elementAt(istart);
			if(jtmp.type==false)
				break;
		}
		if(istart==Crossing.size())
			return "Error generating Dowker notation:\n No undercrossing found.\n";
		for(i=0;i<Crossing.size();i++)
		{
			JCross jtmpi=(JCross)Crossing.elementAt(i);
			JCross jtmpk=(JCross)(jtmpi.partnerObj);
			jtmpi.partner=Crossing.indexOf(jtmpk);
/*
			// Search for partner of crossing i (forms a pair, one odd, one even)
			for(k=0;k<Crossing.size();k++)
			{
				JCross jtmpk=(JCross)Crossing.elementAt(k);
				if(k!=i)
				{
					if(jtmpk.over==jtmpi.over && jtmpk.under==jtmpi.under)
					{
						//partner[i]=k;
						jtmpi.partner=k;
						// assert that jtmpi.type and jtmpk.type should be opposite
						//if(jtmpk.type==jtmpi.type)
						//	return "Error generating Dowker notation:\n Type should be opposite\n";
						break;
					}
				}
			}
			if(k==Crossing.size())
				return "Error generating Dowker notation:\n Incomplete partner found.\n";
*/
		}

		int[] partner2=new int[Crossing.size()];	// have to shift first crossing to first undercrossing
		boolean[] isover=new boolean[Crossing.size()];
		for(i=0;i<Crossing.size();i++)
		{
			JCross jtmpi=(JCross)Crossing.elementAt((istart+i)%Crossing.size());
			//int p=partner[(istart+i)%Crossing.size()];
			int p=jtmpi.partner;
			if(p>=istart)
			{
				p-=istart;
			}
			else
			{
				p-=(istart-Crossing.size());
			}
			partner2[i]=p;
			isover[i]=jtmpi.type;
		}
		String tmp="";

		//for(i=0;i<Crossing.size();i++)
		//{
		//	tmp+=" "+i+" "+partner[i]+"\n";
		//}

		for(i=0;i<Crossing.size();i+=2)	// display even numbered crossings
		{
			int p=partner2[i]+1;	// remember partner,partner2 start at 0
			if(isover[i])
				p=-p;
			tmp+=" "+p+",";
		}

		return "Begin Dowker Representation:\n"+tmp+"\nEnd Dowker Representation\n";
	}	// GenerateDowker


	// Idea based on Vologodskii, et al
	String GenerateAlexander()
	{
		int istart,i,j,k;
		for(istart=0;istart<Crossing.size();istart++)
		{
			JCross jtmp=(JCross)Crossing.elementAt(istart);
			if(jtmp.type==false)
				break;
		}
		if(istart==Crossing.size())
			return "Error generating Alexander matrix:\n No undercrossing found.\n";

		// Determine underpass number and generator number (arc number)
		int undernum=0;
		int gennum=1;
		int numunderpasses=Crossing.size()/2;
		for(i=0;i<Crossing.size();i++)
		{
			JCross jtmpi=(JCross)Crossing.elementAt((istart+i)%Crossing.size());
			if(jtmpi.type==false)
			{
				jtmpi.underpassnum=++undernum;
				if(gennum>=numunderpasses)	// number of arcs equals number of crossings (of one type)
				{
					gennum=1;
				}
				else
				{
					gennum++;
				}
			}
			else	// overpass part of generator or arc
			{
				jtmpi.generatornum=gennum;
			}
		}

		// Write elements of the Alexander matrix into a string
		int writhe=0;
		//int alexmatrix[][]=new int[numunderpasses][numunderpasses];
		BigInteger alexmatrix[][]=new BigInteger[numunderpasses][numunderpasses];
		int tvar;
		try
		{
			tvar=Integer.parseInt(tfAlexVar.getText());
		}
		catch(Exception ex)
		{
			tvar=1;
		}
		String amatrix="";
		for(int l=0;l<Crossing.size();l++)	// twice the number of intersections
		{
			JCross jtmpk=(JCross)Crossing.elementAt((istart+l)%Crossing.size());
			// generate a row of the matrix corresponding to an underpass
			if(jtmpk.type==false)	// only crosstypes for underpasses filled up?
			{
				// Get overpassing generator number of the kth underpass
				//JCross ji=(JCross)Crossing.elementAt(jtmpk.partner);
				JCross ji=(JCross)(jtmpk.partnerObj);
				//if(ji.type!=true) or (jtmpk.type==ji.type) then error
				i=ji.generatornum;
				k=jtmpk.underpassnum;
				if(i==k || i==k+1)	// Rule 1
				{
					for(int m=1;m<=numunderpasses;m++)
					{
						if(m==k) { amatrix+="-1,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(-1); }
						else if(m==k+1) { amatrix+="1,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(1); }
						else { amatrix+="0,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(0); }
					}
					amatrix+="\n";

					// help determine writhe
					// Type II underpass
					if((jtmpk.overdir=='F' && jtmpk.underdir=='L')
						|| (jtmpk.overdir=='R' && jtmpk.underdir=='F')
						|| (jtmpk.overdir=='L' && jtmpk.underdir=='B')
						|| (jtmpk.overdir=='B' && jtmpk.underdir=='R'))
					{
						jtmpk.crosstype=+1;
						writhe+=jtmpk.crosstype;
					}
					// Type I underpass
					if((jtmpk.overdir=='L' && jtmpk.underdir=='F')
						|| (jtmpk.overdir=='F' && jtmpk.underdir=='R')
						|| (jtmpk.overdir=='B' && jtmpk.underdir=='L')
						|| (jtmpk.overdir=='R' && jtmpk.underdir=='B'))
					{
						jtmpk.crosstype=-1;
						writhe+=jtmpk.crosstype;
					}
				}
				else	// Rule 2
				{
					// Type II underpass.  Convention switched 02/11/02, right-handed to left-handed, corresponding to reversing x-direction
					if((jtmpk.overdir=='F' && jtmpk.underdir=='L')
						|| (jtmpk.overdir=='R' && jtmpk.underdir=='F')
						|| (jtmpk.overdir=='L' && jtmpk.underdir=='B')
						|| (jtmpk.overdir=='B' && jtmpk.underdir=='R'))
					{
						jtmpk.crosstype=+1;
						writhe+=jtmpk.crosstype;
						for(int m=1;m<=numunderpasses;m++)
						{
							if(m==k) { amatrix+="-t,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(-tvar); }
							else if(m==k+1)  { amatrix+="1,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(1); }
							else if(m==i) { amatrix+="t-1,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(tvar-1); }
							else { amatrix+="0,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(0); }
						}
						amatrix+="\n";
					}
					// Type I underpass
					if((jtmpk.overdir=='L' && jtmpk.underdir=='F')
						|| (jtmpk.overdir=='F' && jtmpk.underdir=='R')
						|| (jtmpk.overdir=='B' && jtmpk.underdir=='L')
						|| (jtmpk.overdir=='R' && jtmpk.underdir=='B'))
					{
						jtmpk.crosstype=-1;
						writhe+=jtmpk.crosstype;
						for(int m=1;m<=numunderpasses;m++)
						{
							if(m==k) { amatrix+="1,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(1); }
							else if(m==k+1) { amatrix+="-t,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(-tvar); }
							else if(m==i) { amatrix+="t-1,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(tvar-1); }
							else { amatrix+="0,"; alexmatrix[k-1][m-1]=BigInteger.valueOf(0); }
						}
						amatrix+="\n";
					}
				}
			}
		}

		String tmp="Begin Alexander Matrix\n"+amatrix+"End Alexander Matrix\nWrithing Number: "+writhe+"\n";

		// Calculate Determinant of n-1 minor of the Alexander Matrix.  This is almost the Alexander Polynomial evaluated at t.
		if(cbDet.getState())
		{
			if(numunderpasses>1)
				numunderpasses--;	// n-1 minor
/*** test
			alexmatrix=new BigInteger[3][3];
			numunderpasses=3;
			alexmatrix[0][0]=BigInteger.valueOf(-14);alexmatrix[0][1]=BigInteger.valueOf(4);alexmatrix[0][2]=BigInteger.valueOf(40);
			alexmatrix[1][0]=BigInteger.valueOf(5);alexmatrix[1][1]=BigInteger.valueOf(-2);alexmatrix[1][2]=BigInteger.valueOf(1);
			alexmatrix[2][0]=BigInteger.valueOf(-3);alexmatrix[2][1]=BigInteger.valueOf(6);alexmatrix[2][2]=BigInteger.valueOf(37);
***/
			BigInteger bipivot,bi2,biexcess=BigInteger.valueOf(1);
			for(i=0;i<numunderpasses;i++)
			{
				// get 'pivot' row, first element is nonzero
				for(j=i;j<numunderpasses;j++)
				{
					if(!alexmatrix[j][i].equals(BigInteger.valueOf(0)))
					{
						break;
					}
				}
				// if no pivot row, determinant is zero
				if(j==numunderpasses)
				{
					tmp+="Determinant of n-1 minor is zero\n";
					break;
				}
				// swap
				if(j!=i)
				{
					for(k=i;k<numunderpasses;k++)
					{
						BigInteger tmpbi=new BigInteger(alexmatrix[i][k].toString());
						alexmatrix[i][k]=alexmatrix[j][k];
						alexmatrix[j][k]=tmpbi;
					}
				}
				// subtract off multiple of pivot row from rest of rows
				for(j=i+1;j<numunderpasses;j++)	// for each row
				{
					bipivot=alexmatrix[i][i];
					bi2=alexmatrix[j][i];
					alexmatrix[j][i]=BigInteger.valueOf(0);
					biexcess=biexcess.multiply(bipivot);
					for(k=i+1;k<numunderpasses;k++)	// for each column, perform subtraction (elimination)
					{
						alexmatrix[j][k]=alexmatrix[j][k].multiply(bipivot);
						alexmatrix[j][k]=alexmatrix[j][k].subtract(alexmatrix[i][k].multiply(bi2));
					}
				}
			}
			// take product of diagonal elements and divide by correction factor
			for(i=1;i<numunderpasses;i++)
			{
				alexmatrix[0][0]=alexmatrix[0][0].multiply(alexmatrix[i][i]);
			}
			alexmatrix[0][0]=alexmatrix[0][0].divide(biexcess);
			tmp+="Determinant of n-1 minor (t="+tvar+"): "+alexmatrix[0][0].toString()+"\n";
		}

		return tmp;
	}	// GenerateAlexander

	// Assume p1 and p2 differ in only one coordinate, by 1
	char SegDir(int x1, int y1, int z1, int x2, int y2, int z2)
	{
		if(x1-x2==1) return 'B';	//relative to p1
		if(x2-x1==1) return 'F';
		if(y1-y2==1) return 'L';
		if(y2-y1==1) return 'R';
		if(z1-z2==1) return 'D';
		if(z2-z1==1) return 'U';

		return 'N';	//error?
	}

	// Emulate ParsePolymer.  Gauss code for knot as input.
	void ParseGauss(String code)	//This feature came out from C++ code first
	{
		Crossing.removeAllElements();	// clear Crossing first

		StringTokenizer st=new StringTokenizer(code,",\r\n");
		int i=0;
		while(st.hasMoreTokens())
		{
			i++;
			String sc=st.nextToken();
			sc=sc.trim();
			if(sc.length()<3)
			{
				taOutput.append("Error parsing Gauss code:\n Incorrect input at element "+i+"\n");
				return;
			}
			boolean ab;
			if(sc.charAt(0)=='a')
				ab=true;
			else if(sc.charAt(0)=='b')
				ab=false;
			else
			{
				taOutput.append("Error parsing Gauss code:\n Incorrect a-b input at element "+i+"\n");
				return;
			}
			int num;
			try
			{
				num=Integer.parseInt(sc.substring(2));
			}
			catch(Exception ex)
			{
				taOutput.append("Error parsing Gauss code:\n Incorrect input label at element "+i+"\n");
				return;
			}
			char odir='L',udir='F';
			if(sc.charAt(1)=='-')
				num=-num;
			if(num>0)
			{
				odir='F';
				udir='L';
			}
			//Crossing.push_back(new CCross(ab,pc->num,pc->num,0,0,true,odir,udir));
			Crossing.addElement(new JCross(ab,num,num,-1,-1,true,odir,udir));
		}

		if(Crossing.size()>1)
		{
			taOutput.append(MatchPartnerObjects());
			taOutput.append(GenerateDowker());
			taOutput.append(GenerateAlexander());
		}
	}

	/**
	 * The entry point for the applet. 
	 */
	public void init()
	{
		initForm();
		usePageParams();
	}

	/**
	 * Reads parameters from the applet's HTML host and sets applet
	 * properties.
	 */
	private void usePageParams()
	{
	}

	/**
	 * Intializes values for the applet and its components
	 */
	void initForm()
	{
		this.setLayout(new GridLayout(2,2,5,5));
		this.setBackground(Color.black);

		cv3D=new J3D(this.getSize().width/2);
		cv3D.polyref=Polymer;
		cv3D.setBackground(Color.black);
		this.add(cv3D);

		cv2D.setBackground(Color.blue);
		cv2D.polyref=Polymer;
		cv2D.crossref=Crossing;
		this.add("Center",cv2D);	//Center?

		Panel p,p2,p3;
		p=new Panel();
		p.setLayout(new BorderLayout());
		p.setBackground(Color.gray);
		p2=new Panel();	//panel for control buttons
		p2.add(new Label("3D Perspective"));
		p2.add(btUp);
		btUp.addActionListener(new BL(btUp.getLabel()));
		p2.add(btDown);
		btDown.addActionListener(new BL(btDown.getLabel()));
		p2.add(btLeft);
		btLeft.addActionListener(new BL(btLeft.getLabel()));
		p2.add(btRight);
		btRight.addActionListener(new BL(btRight.getLabel()));
		p3=new Panel();
		p3.setBackground(Color.black);
		p3.add(btInput);
		chInput=new Choice();
		chInput.add("x y z");
		chInput.add("z x y");
		chInput.add("y z x");
		chInput.add("Gauss");
		p3.add(chInput);
		p3.add(btClearInput);
		btClearInput.addActionListener(new BL(btClearInput.getLabel()));
		btInput.addActionListener(new BL(btInput.getLabel()));
		p.add("North",p2);
		p.add("Center",p3);
		p.add("South",taInput);
		this.add(p);
// Initial data
taInput.setText("0 0 0\n0 0 1\n0 0 2\n0 0 3\n1 0 3\n1 0 2\n1 0 1\n1 0 0\n2 0 0\n2 0 1\n2 0 2\n2 0 3\n2 1 3\n2 1 2\n2 2 2\n2 2 1\n2 2 0\n2 3 0\n3 3 0\n3 3 1\n2 3 1\n2 3 2\n3 3 2\n3 3 3\n2 3 3\n1 3 3\n1 3 2\n1 3 1\n1 3 0\n0 3 0\n0 3 1\n0 3 2\n0 3 3\n0 2 3\n0 1 3\n1 1 3\n1 1 2\n1 1 1\n1 1 0\n2 1 0\n2 1 1\n3 1 1\n3 1 2\n3 1 3\n3 0 3\n3 0 2\n3 0 1\n3 0 0\n3 1 0\n3 2 0\n3 2 1\n3 2 2\n3 2 3\n2 2 3\n1 2 3\n1 2 2\n1 2 1\n1 2 0\n0 2 0\n0 2 1\n0 2 2\n0 1 2\n0 1 1\n0 1 0\n");


//taInput.setText("0 0 0\n0 0 1\n0 1 1\n0 1 0\n0 2 0\n0 3 0\n1 3 0\n2 3 0\n2 3 1\n1 3 1\n1 3 2\n1 3 3\n2 3 3\n2 3 2\n3 3 2\n3 3 3\n3 2 3\n3 2 2\n2 2 2\n2 2 3\n2 1 3\n3 1 3\n3 0 3\n2 0 3\n1 0 3\n0 0 3\n0 1 3\n0 1 2\n0 0 2\n1 0 2\n1 1 2\n1 1 3\n1 2 3\n1 2 2\n0 2 2\n0 2 3\n0 3 3\n0 3 2\n0 3 1\n0 2 1\n1 2 1\n1 2 0\n2 2 0\n3 2 0\n3 3 0\n3 3 1\n3 2 1\n2 2 1\n2 1 1\n3 1 1\n3 1 2\n2 1 2\n2 0 2\n3 0 2\n3 0 1\n2 0 1\n2 0 0\n3 0 0\n3 1 0\n2 1 0\n1 1 0\n1 1 1\n1 0 1\n1 0 0\n");
//taInput.setText("0 0 0\n0 0 1\n1 0 1\n1 0 0\n1 1 0\n1 1 1\n0 1 1\n0 1 0\n");

		p=new Panel();
		p.setLayout(new BorderLayout());
		p.setBackground(Color.gray);
		p2=new Panel();
		p2.setBackground(Color.blue);
		p2.add(cbReduce);
		p2.add(cbDet);
		tfAlexVar.setBackground(Color.white);
		p2.add(tfAlexVar);
		p2.add(btClear);
		btClear.addActionListener(new BL(btClear.getLabel()));
		p3=new Panel();
		p3.add(new Label("Plane Projection, x -> x + z/side, y -> y + z/side"));
		p.add("Center",p2);
		p.add("North",p3);
		p.add("South",taOutput);
		this.add(p);
	}

	class BL implements ActionListener
	{
		String ButtonLabel;

		BL(String s)
		{
			ButtonLabel=s;
		}

		public void actionPerformed(ActionEvent e)
		{
			if(ButtonLabel.equals("Process Input"))
			{
				int c=chInput.getSelectedIndex();
				if(c==3)
				{
					Polymer.size=-1;	//help prevent drawing
					//cv3D.Reset();
					//cv3D.repaint();
					taOutput.append("Processing Gauss code...\n");
					ParseGauss(taInput.getText());
					//cv2D.repaint();
					taOutput.append("...Done\n");
				}
				//else if(c==0,1,2)
				else
				{
					if(Polymer.Parse(taInput.getText(),c))
					{
						taOutput.append("Processing Coordinates...\n");
						//cv3D.Reset();
						//cv3D.repaint();
						//cv2D.repaint();
						try{ ParsePolymer(); }catch(Exception e2){taOutput.append("Exception:"+e2+"\n");}
						taOutput.append("...Done\n");
					}
					else
					{
						Polymer.size=0;
						//cv3D.Reset();
						//cv3D.repaint();
						//cv2D.repaint();
						taOutput.append("Error! Please check your input.\n");
					}
				}
				cv3D.Reset();
				cv3D.repaint();
				cv2D.repaint();
			}
			else if(ButtonLabel.equals("Clear Output"))
			{
				taOutput.setText("");
			}
			else if(ButtonLabel.equals("Clear Input"))
			{
				taInput.setText("");
			}
			else
			{
				cv3D.HandleEvent(ButtonLabel);
				//cv2D.repaint();
			}
		}
	}

	// Need this if compiled as standalone app
	public static void main(String[] args)
	{
		JKnot ja=new JKnot();
		ja.setSize(400,400);
		ja.init();
		JFrame jf=new JFrame("Knots");
		jf.setLayout(new GridLayout(1,1));
		jf.setSize(410,410);
		jf.add(ja);
		jf.show();
	}
}

// Need this if compiled as standalone app
class JFrame extends Frame implements WindowListener
{
	public JFrame()
	{
		super();
		addWindowListener( this );
	};

	public JFrame(String title)
	{
		super(title);
		addWindowListener( this );
	};

	public void windowClosing( WindowEvent evt )
	{
		((JKnot)getComponent(0)).destroy();
		System.exit(0);
	};

	public void windowOpened( WindowEvent evt )
	{
	};

	public void windowClosed( WindowEvent evt )
	{
	};

	public void windowDeiconified( WindowEvent evt )
	{
	};

	public void windowActivated( WindowEvent evt )
	{
	};

	public void windowIconified( WindowEvent evt )
	{
	};

	public void windowDeactivated( WindowEvent evt )
	{
	};
}
