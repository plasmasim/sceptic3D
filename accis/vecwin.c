#include <windows.h>
#include <stdio.h>

#define accis_path_max 4000
/*  typedef struct {  */
/*    int x; */
/*    int y; */
/*  } XPoint; */
/*  XPoint accis_path[accis_path_max]; */
  POINT accis_path[accis_path_max];
int accis_pathlen=0;

#define accis_maxPixels 16
COLORREF accis_pixels[accis_maxPixels]=
{ 0x00FFFFFF,
  0x00770000,  
  0x00007700,  
  0x00777700,  
  0x000000AA,  
  0x00770077,  
  0x000077AA,
  0x00777777,  
  0x00AAAAAA,  
  0x00FF0000,  
  0x0000FF00,  
  0x00FFFF00,  
  0x000000FF,  
  0x00FF00FF,  
  0x0000FFFF,  
  0x00000000,  
};

/* Unused on windows.
char *accis_colornames[accis_maxPixels]=
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
*/

#define ACCIS_DEFWIDTH 640
#define ACCIS_DEFHEIGHT 480

static char g_szClassName[] = "MyWindowClass";
static HINSTANCE g_hInst = NULL;
static int wWidth, wHeight; /* Window width and height */
static RECT wRect;
PAINTSTRUCT ps;
HDC hdcMemory, hdcWindow;
BITMAP bm;
HBITMAP hbm; 
int first=1;
static HBRUSH hbr;
static HBRUSH accis_brush;
static HPEN accis_pen;
  
static WNDCLASSEX WndClass;
static HWND hwnd;
static MSG Msg;

LRESULT CALLBACK WndProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{

   switch(Message){
   case WM_CREATE:
     break;
   case WM_SIZE:  
     /* This is called with 0,0 when iconized. That is bad. So don't resize
	now.*/
     wWidth = LOWORD(lParam);  // width of client area 
     wHeight = HIWORD(lParam); // height of client area 
/*       printf("Resized: %d %d\n",wWidth,wHeight); */
     /*     AccisRenewBitmap(hwnd);*/
     break;
   case WM_PAINT:
/*       printf("Paint BitBlt\n"); */
     /*Trying for private DC*/
/*       hdcWindow = BeginPaint(hwnd, &ps); */
     BitBlt(hdcWindow, 0, 0, wWidth, wHeight, hdcMemory, 0, 0, SRCCOPY);
/*         EndPaint(hwnd, &ps); */
     GetClientRect(hwnd,&wRect);
     ValidateRect(hwnd,&wRect);/* Needed without endpaint. Not perfect.
				Resize leads to continuous BitBlt.
			       Needs new/current window size.*/
     break;
   case WM_CLOSE:
       /* It is not safe to do the following from the right-hand X.
	  Makes the program hang. The Default is also bad so act the same.*/
/*       DestroyWindow(hwnd); */
/*       break; */
   case WM_DESTROY:
   case WM_LBUTTONDOWN: /* Click in window to exit from this message loop. */
     DeleteDC(hdcMemory);
     PostQuitMessage(0);
     break;
   default:
     return DefWindowProc(hwnd, Message, wParam, lParam);
   }
   return 0;
}


/**********************************************************************/
/* Delete the old bitmap and make a new one. */
void AccisRenewBitmap(HWND hwnd){
/*      WPARAM wParam; */
/*      LPARAM lParam; */
     if(!first){
/*         printf("Deleting hdcMemory and hbm\n"); */
       DeleteDC(hdcMemory);
       DeleteObject(hbm);
     }else{
       first=0;
     }
     GetClientRect(hwnd, &wRect);
     hdcWindow = GetDC(hwnd);
     /* Hoped this is not nec when we use GetDC instead. No worse.*/
/*         hdcWindow = BeginPaint(hwnd, &ps);   */
     hdcMemory = CreateCompatibleDC(hdcWindow);            
/*       printf("Creating bitmap, %d %d\n",wWidth,wHeight); */
     /* This needs to use the window device context to get color.*/
     hbm=CreateCompatibleBitmap(hdcWindow,wWidth,wHeight);
/*         EndPaint(hwnd, &ps);  */
     SelectObject(hdcMemory,hbm);
     GetObject(hbm,sizeof(bm),&bm);
     /* Have to set bitmap background. Not set automatically.*/
     hbr=CreateSolidBrush(0x00FFFFFF);
     /* Removed this to rely on WM_PAINT which ought to be called auto.*/
/*       BitBlt(hdcWindow, 0, 0, wWidth, wHeight,hdcMemory,0,0,SRCCOPY); */
     FillRect(hdcMemory,&wRect,hbr);
     InvalidateRect(hwnd,&wRect,TRUE); /* Force whole window draw */
     /*This does not seem to be necessary once we invalidate */
/*         WndProc(hwnd,WM_PAINT,wParam,lParam);  */
}

/* svga does the job of WinMain setting up the window etc.*/

int svga_(int *scrxpix, int *scrypix, int *vmode, int *ncolor)
{
  static int second=0;
  /*  extern int f__xargc;
      extern char **f__xargv;
      int svga_argc=0;
      char *svga_argv[1]; 
  */
  /* For some bizarre reason the following two definitions are important
     even though they do not seem to be used. Without them the window does
     not pop up.
  */
  HINSTANCE hInstance;
  HINSTANCE hPrevInstance; 
  LPSTR lpCmdLine; 
  /*nCmdShow tells how to open window */
  int nCmdShow=SW_SHOWNORMAL;

  g_hInst=GetModuleHandle(NULL);
/*    printf("Entered svga\n"); */
  if(second == 0){
    WndClass.cbSize        = sizeof(WNDCLASSEX);
    WndClass.style         = CS_OWNDC;
    WndClass.lpfnWndProc   = WndProc;
    WndClass.cbClsExtra    = 0;
    WndClass.cbWndExtra    = 0;
    WndClass.hInstance     = g_hInst;
    WndClass.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    WndClass.hCursor       = LoadCursor(NULL, IDC_ARROW);
    WndClass.hbrBackground = (HBRUSH)(COLOR_BTNFACE+1);
    WndClass.lpszMenuName  = NULL;
    WndClass.lpszClassName = g_szClassName;
    WndClass.hIconSm       = LoadIcon(NULL, IDI_APPLICATION);
    
    if(!RegisterClassEx(&WndClass))
      {
	MessageBox(0, "Window Registration Failed!", "Error!",
		   MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
	return 0;
      }
    
    hwnd = CreateWindowEx(
			  WS_EX_CLIENTEDGE,
			  g_szClassName,
			  "Accis Window",
			  WS_OVERLAPPEDWINDOW,
			  CW_USEDEFAULT, CW_USEDEFAULT,
			  ACCIS_DEFWIDTH+12, ACCIS_DEFHEIGHT+31,
			  NULL, NULL, g_hInst, NULL);
    
    if(hwnd == NULL)
      {
	MessageBox(0, "Accis Window Creation Failed!", "Error!",
		   MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
	return 0;
      }
    
/*      printf("About to show window\n"); */
    ShowWindow(hwnd, nCmdShow);
/*      printf("About to update window\n"); */
    UpdateWindow(hwnd);

    second++;
  }

  /* Originally this was implicit in resize. mingw did not like it.*/
  AccisRenewBitmap(hwnd);
  
/*    printf("Entering peek message loop. Second=%d\n",second); */
  while(PeekMessage(&Msg, hwnd, 0, 0, PM_REMOVE))
    {
      TranslateMessage(&Msg);
      DispatchMessage(&Msg);
    }
  hdcWindow = GetDC(hwnd); /*Ensure we have still selected the window.*/

  *vmode=88;
  *ncolor=15;
  *scrxpix=wWidth;
  *scrypix=wHeight;
  
  return Msg.wParam;
}
/*************************************************************************/
int txtmode_()
{
  while(GetMessage(&Msg, hwnd, 0, 0))
    {
      TranslateMessage(&Msg);
      DispatchMessage(&Msg);
    }
  return 0;
}
/* ******************************************************************** */
int vec_(long *px,long *py,long *ud)
{ /*  Draw vector on screen, with pen up or down. */
    static int px1=0,py1=0,px2=0,py2=0;
    extern POINT accis_path[];
    extern int accis_pathlen;

    px1=px2;
    py1=py2;
    px2 = *px;
    py2 = *py;
    if( *ud != 0) {
/*        MoveToEx(hdcWindow,px1,py1,NULL); */
/*        MoveToEx(hdcMemory,px1,py1,NULL); */
      LineTo(hdcWindow,px2,py2);
      LineTo(hdcMemory,px2,py2);
      if(accis_pathlen<accis_path_max){      /* Add point to path */
	accis_pathlen++;
      }
    }else{ /* Restart path */
      accis_pathlen=0;
      MoveToEx(hdcWindow,px2,py2,NULL);
      MoveToEx(hdcMemory,px2,py2,NULL);
    }
    accis_path[accis_pathlen].x=*px;
    accis_path[accis_pathlen].y=*py;

    return 0;
} /* vec_ */
/* ******************************************************************** */

int scolor_(long *li)
{
  if((*li < accis_maxPixels) && (*li >= 0)){
    DeleteObject(accis_pen);
    accis_pen=CreatePen(PS_SOLID,0,accis_pixels[(int)(*li)]);
    DeleteObject(accis_brush);
    accis_brush=CreateSolidBrush(accis_pixels[(int)(*li)]);
    SelectObject(hdcWindow,accis_pen);
    SelectObject(hdcMemory,accis_pen);
    SelectObject(hdcWindow,accis_brush);
    SelectObject(hdcMemory,accis_brush);
    return 1;
  }else{    
    return 0;
  }
}
/* ********************************************************************* */
int vecfill_()
{
  SetPolyFillMode(hdcWindow,WINDING);
  SetPolyFillMode(hdcMemory,WINDING);
  Polygon(hdcWindow,accis_path,accis_pathlen+1);
  Polygon(hdcMemory,accis_path,accis_pathlen+1);
  return 0;
}
/* ********************************************************************* 
int main()
{
  long thecolor=0;
  int scrxpix, scrypix, vmode, ncolor;
  long x,y,pen;
  int i,j;
  svga_(&scrxpix,&scrypix, &vmode, &ncolor);
  printf("Values returned from svga:\n %d %d %d %d\n", scrxpix, scrypix, vmode, ncolor); 
  x=0;y=100;
  for (i=1;i<16;i++){
    y=y+5;
    thecolor=i;
    scolor_(&thecolor);
    for (j=1;j<6;j++){
      pen=0;
      x=0;
      vec_(&x,&y,&pen);
      x=100;pen=1;
      vec_(&x,&y,&pen);
      y++;
    }
  }
  txtmode_();
  return 0;
}
*/
/* Notes.
   It is something to do with beginpaint, endpaint being omitted. Looks
   as if one can't do that. Seems as if you have to do beginpaint anew
   every time otherwise it paints in the wrong place.

   11 Jul 2001
   The problem with missing parts of the plot occurs if we do a click in
   window for a continuation when the window is partly hidden. This must
   be a timing problem. There must be some kind of thread that is not 
   completing before the bitblt happens.

   This now seems to be fixed by putting all blts into the WM_PAINT and
   making sure it gets called by the invalidaterect.

   It seems that we need to study Display Device Contexts and figure out how
   to use a Privaate Device Context that retains its graphic objects.
 */

/* Rotation code*/
float xeye,yeye,zeye;
float xeye0,yeye0,zeye0;
float accis_x0,accis_y0;

void extern butdown_();
void extern viewrot_();
void extern cubeupd_();
void extern butup_();

void accis_butdown(lMsg)
     MSG lMsg;
{ /* accis_butdown */
  accis_x0=LOWORD(lMsg.lParam);
  accis_y0=HIWORD(lMsg.lParam);
  butdown_(&xeye0,&yeye0,&zeye0);
  xeye=xeye0; yeye=yeye0; zeye=zeye0;
  }

int cmpmessage(UINT message)
{
   GetMessage(&Msg, hwnd, 0, 0);
   if(Msg.message == message){   return 1;   }else{  return 0;   }
}

int eye3d_(value)
     int *value;
{
  float xmoved,ymoved;
    while(!cmpmessage(WM_LBUTTONDOWN )){
      TranslateMessage(&Msg);
      DispatchMessage(&Msg);
      }
/*      printf("Got first button press...");  */
  accis_butdown(Msg);
  while(!cmpmessage(WM_LBUTTONUP)){
      /* If we started moving without setting eye nonzero, set it first*/
      if(xeye0==0. && yeye0==0. && zeye0==0.) accis_butdown(Msg); 
      xmoved=LOWORD(Msg.lParam)-accis_x0;
      ymoved=HIWORD(Msg.lParam)-accis_y0;
/*          printf("eye: %f %f %f\n",xeye,yeye,zeye);  */
      viewrot_(&xmoved,&ymoved,&xeye0,&yeye0,&zeye0,&xeye,&yeye,&zeye);
      cubeupd_(&xeye,&yeye,&zeye);
      {/* Full redraw section not working.
	PostMessage(hwnd,WM_LBUTTONDOWN,0,0);
	break; */
      }
  }
/*      printf("Exited eye3d loop\n");  */
  butup_(&xeye,&yeye,&zeye);
  *value=0;
  if(  !( xeye==xeye0 && yeye==yeye0 && zeye==zeye0) ) *value=1;
  return 0;
}
/* 10 May 2003 got rotation working.*/
