/* Xwindow driver for accis plotting */
/* Fortran callable (f2c) routines */
/* ********************************************************************** */
/*
Refreshing version.
*/
  #include  <X11/StringDefs.h> 
  #include  <X11/Intrinsic.h> 
  #include  <X11/Core.h> 

/* Globals for these routines*/
Widget accis_wshell;
Widget accis_drawing;
Display *accis_display;
Drawable accis_window;
Pixmap accis_pixmap;

GC accis_gc;
int accis_depth;
Colormap accis_colormap;

struct Screen_Size {  
    Dimension width; /* unsigned short */
    Dimension height;
  };
    static struct Screen_Size s_s;
/* Static maximum number of points in path */
#define accis_path_max 4000
XPoint accis_path[accis_path_max];
int accis_pathlen=0;

#define a_maxPixels 16
unsigned long accis_pixels[a_maxPixels];

char *accis_colornames[a_maxPixels]=
{
  "White",
  "MediumBlue",
  "SeaGreen",
  "MediumTurquoise",
  "Firebrick",
  "Orchid",
  "Brown",
  "LightGrey",
  "SlateGrey",
  "Blue",
  "Green",
  "Cyan",
  "Red",
  "Magenta",
  "Yellow",
  "Black"
};



/* Subroutine */ 
int svga_(scrxpix, scrypix, vmode, ncolor)
int *scrxpix, *scrypix, *vmode, *ncolor;
{  
  static int second=0;
  extern int xargc;
  extern char **xargv;
  int svga_argc=0;
  char *svga_argv[1]; 
  static int n;
  static Arg wargs[10];
  int theDepth;
  Colormap theColormap;
  extern void config_handler();
  extern void accis_refresh();

  if(second == 0){
    accis_wshell = XtInitialize("accis","Accis", NULL, 0,
      &xargc, xargv);  /* This refers to the f2c command line args.
			It may only work with f2c, therefore. */
      /* &svga_argc, svga_argv); alternate */

    accis_drawing = XtCreateManagedWidget("drawing",coreWidgetClass,
		 accis_wshell, NULL, 0);

    *vmode=88;
    *ncolor=15;
    /* Set up a default size of the drawing window widget. 
       This is overruled by Accis*geometry resources setting wshell size.
       */
    n = 0;
    XtSetArg(wargs[n], XtNheight, 480); n++;
    XtSetArg(wargs[n], XtNwidth, 640); n++;
    XtSetValues(accis_drawing, wargs, n);
    XtRealizeWidget(accis_wshell);
    accis_display = XtDisplay(accis_drawing);
    accis_window = XtWindow(accis_drawing);
    accis_gc = XCreateGC(accis_display, accis_window, 0, NULL);
    accis_depth=DefaultDepth(accis_display,0);
    accis_colormap=DefaultColormap(accis_display,0);
    initDefaultColors();
    /* Leave setup for resizing. */
    n = 0;
    XtSetArg(wargs[n], XtNheight, &s_s.height); n++;
    XtSetArg(wargs[n], XtNwidth, &s_s.width); n++;
    /* Pixmap setup */
    XtGetValues(accis_wshell, wargs, n);
    accis_pixmap=XCreatePixmap(accis_display,accis_window,
			       s_s.width,s_s.height,accis_depth);
    XtAddEventHandler(accis_drawing,ExposureMask,FALSE,
		      accis_refresh,NULL);
    XSelectInput(accis_display,accis_window,
		 KeyPressMask | ExposureMask | ButtonPress
		 | FocusChangeMask | EnterWindowMask ); /* Tell events */
    second++;
  }else{
    /* Need to have this to make the get correct. */
    ManageEvents();
  } 
  /* This doesn't need a configuration event handler to work correctly. */
  XtGetValues(accis_wshell, wargs, n);
  *scrxpix=s_s.width;
  *scrypix=s_s.height;
  XFlush(accis_display);
  XClearWindow(accis_display,accis_window);
  XSetForeground(accis_display,accis_gc,accis_pixels[0]);  
  /* Clear the window to background color. VMS doesn't do correctly.
     But also this seems to fix the bad match error. */
  XFillRectangle(accis_display,accis_window,accis_gc,0,0,
		 s_s.width,s_s.height);
  /* Clear the pixmap  */
  XFillRectangle(accis_display,accis_pixmap,accis_gc,0,0,
		 s_s.width,s_s.height);
  XSetForeground(accis_display,accis_gc,accis_pixels[15]); 
  XSetInputFocus(accis_display, accis_window, RevertToPointerRoot,
      CurrentTime);
  return 0;
}

/* ******************************************************************** */
initDefaultColors()
{
  XColor theRGBColor,theHardColor;
  int status;
  int i;
  for(i=0;i<a_maxPixels;i++){
    status=XLookupColor(accis_display,accis_colormap,accis_colornames[i],
			&theRGBColor,&theHardColor);
    if(status !=0){
      status=XAllocColor(accis_display,accis_colormap,&theHardColor);
      if(status !=0){
	accis_pixels[i]=theHardColor.pixel;
      }else{
	accis_pixels[i]=BlackPixel(accis_display,0);
      }
    }else{
      accis_pixels[i]=BlackPixel(accis_display,0);
    }
  }
}

/* ******************************************************************** */
/* End plotting and return to text editing. */
/* Subroutine */ 
/* #include <curses.h>*/
int txtmode_()
{
  XEvent event; 
  XFlush(accis_display);
  do{
    /*    printf("Executing XtNextEvent"); */
    XtNextEvent(&event);
    /* XNextEvent(accis_display,&event); is equivalent */
    XtDispatchEvent(&event);  
    /*    printf("The event type: %d\n",event); */
  }while(event.type != ButtonPress && event.type != KeyPress );
  /* Here we should give the focus back to parent, but I don't see how.
    XSetInputFocus(accis_display, ??, PointerRoot,
		   CurrentTime); */
}

/* ********************************************************************* */
/* Subroutine */ int scolor_(li)
long *li;
{
  /* *ncolor=*li; */
  if((*li < a_maxPixels) && (*li >= 0)){
    XSetForeground(accis_display,accis_gc,accis_pixels[(int) *li]);
    return 1;
  }else{    
    return 0;
  }
} /* scolor_ */

/* ******************************************************************** */
/* Subroutine */ int vec_(px, py, ud)
long *px, *py, *ud;
{ /*  Draw vector on screen, with pen up or down. */
    static int px1=0,py1=0,px2=0,py2=0;
    extern XPoint accis_path[];
    extern int accis_pathlen;

    px1=px2;
    py1=py2;
    px2 = *px;
    py2 = *py;
    if( *ud != 0) {
      XDrawLine(XtDisplay(accis_drawing),XtWindow(accis_drawing), accis_gc,
		  px1,py1,px2,py2);
      XDrawLine(XtDisplay(accis_drawing),accis_pixmap, accis_gc,
		  px1,py1,px2,py2);
      if(accis_pathlen<accis_path_max){      /* Add point to path */
	accis_pathlen++;
      }
    }else{ /* Restart path */
      accis_pathlen=0;
    }
    accis_path[accis_pathlen].x=*px;
    accis_path[accis_pathlen].y=*py;
/*    XFlush(accis_display);
 Flush removed here. Now relies on txtmode to flush display.  */
    return 0;
} /* vec_ */
/* ******************************************************************** */
int vecfill_()
{
    extern XPoint accis_path[];
    extern int accis_pathlen;
    if(accis_pathlen>1){ /* If path is more than 2 points, fill. */
      XFillPolygon(accis_display,accis_window,accis_gc,
		   accis_path,accis_pathlen+1,Nonconvex,CoordModeOrigin);
      XFillPolygon(accis_display,accis_pixmap,accis_gc,
		   accis_path,accis_pathlen+1,Nonconvex,CoordModeOrigin);
    }
}
/* ******************************************************************** */
/* Not currently in use. */
void config_handler(w,cs_s,event)
Widget w;
struct Screen_Size cs_s;
XEvent *event;
{
cs_s.width=event->xconfigure.width;
cs_s.height=event->xconfigure.height;
printf("config_handler width= %d, height= %d\n",cs_s.width,cs_s.height);
}
/* ******************************************************************** */
/* ******************************************************************** */
void accis_refresh(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  XCopyArea(XtDisplay(w),accis_pixmap,accis_window,accis_gc,0,0,
	    s_s.width,s_s.height,0,0);
  XFlush(accis_display);
}
/* ******************************************************************** */

ManageEvents()
{
  XEvent event; 
  while(XtPending()){
    XtNextEvent(&event);
    XtDispatchEvent(&event);  
  }
}
/* ******************************************************************** */
/* Testing only */
/*
main()
 {
   long li=1;
   int x=200,y=150,pen=1;
   int a,b,c,d;
   int ch;
   while(ch != 'q'){
     printf("Calling svga\n");
     svga_(&a,&b,&c,&d);
     printf("Returned, x= %d, y= %d, mode= %d, ncolor= %d,\n",
	    a,b,c,d);
     printf("Drawing stuff directly\n");
     draw_graphics(accis_drawing);
     printf("Drawing line using vec_\n");
     for(a=0;a<16;a++){
       li=a;
       scolor_(&li);
       for(b=0;b<10;b++){
	 pen=0;x=200;y=150+10*a+b;
	 vec_(&x,&y,&pen); 
	 pen=1;x=400;y=150+10*a+b;
	 vec_(&x,&y,&pen); 
       }
     }
     ch=getchar();
   }
 }

  draw_graphics(w)
  Widget w; {
       Display *display;
       Drawable window;
       GC gc;
       int store;
       XWindowAttributes attributes;
       XSetWindowAttributes sattributes;
       unsigned long valuemask;

       display = XtDisplay(w);
       window = XtWindow(w);
       store= DoesBackingStore(DefaultScreenOfDisplay(display));
       printf("DoesBackingstore return, %d of 
               NotUseful,WhenMapped,Always %d,%d,%d\n",
	      store,NotUseful,WhenMapped,Always);
       XGetWindowAttributes(display,window,&attributes);
       printf("Got attributes. backing_store,planes,pixel,%d %x %u\n"
	      ,attributes.backing_store,attributes.backing_planes,
	      attributes.backing_pixel);


       gc = XCreateGC(display, window, 0, NULL);
       XSetForeground(display, gc, 50);
       XSetBackground(display, gc, 0);

       XDrawLine(display, window, gc, 10, 10, 400, 400);
       XDrawRectangle(display, window, gc, 75, 110, 150, 100);
       XDrawArc(display, window, gc, 75, 110, 150, 100, 45*64, 120*64);

       XFreeGC(display, gc);
  }


*/

/* ******************************************************************** */

/* This broke with a complaint about pixmaps.

       
       valuemask=CWBackingStore || CWBackingPlanes || CWBackingPixel;
       sattributes.backing_store=2;
       sattributes.backing_planes=255;
       sattributes.backing_pixel=1;
       XChangeWindowAttributes(display, window, valuemask, &sattributes);
       XGetWindowAttributes(display,window,&attributes);
       printf("Got attributes. backing_store,planes,pixel,%d %x %u\n"
	      ,attributes.backing_store,attributes.backing_planes,
	      attributes.backing_pixel);
	      
	      */
