/*
 * Serial implememntation of Mandelbrot program.
 *
 * This program computes and displays the Mandelbrot set.
 * By default, it examines all points in the complex plane
 * that have both real and imaginary parts between -2 and 2.  
 * 
 * Code originally obtained from Web site for Wilkinson and
 * Allen's text on parallel programming:
 *   http://www.cs.uncc.edu/~abw/parallel/par_prog/
 * 
 * Reformatted and revised by B.Massingill.
 * Further reformatted by MEH.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <mpi.h>
//#include "mandelbrot-gui.h"     /* has setup(), interact() */


/* Default values */
#define N           2           /* size of problem space (x, y from -N to N) */
#define NPIXELS     1000        /* size of display window in pixels */


/* Structure definition for complex numbers */
typedef struct {
  double real, imag;
} complex;

/* Structure definition for the Xwindow */
typedef struct {
    Display *display;
    Window win;
    GC gc;
} XDATA;

/* Shorthand for some commonly-used types */
typedef unsigned int uint;
typedef unsigned long ulong;


/* Initialize for graphical display. 
 *   -- Width, height are dimensions of display, in pixels. */
int setup( uint width, uint height, Display **display, Window *win,
           GC *gc, ulong *min_color, ulong *max_color )
{
  /* Variables for graphical display */
  uint x = 0, y = 0;                  /* window position */
  uint border_width = 4;              /* border width in pixels */
  uint disp_width, disp_height;       /* size of screen */
  uint screen;                        /* which screen */

  char *window_name = "Mandelbrot Set", *disp_name = NULL;
  ulong valuemask = 0;
  XGCValues values;

  ulong white, black;                 /* white, black pixel values */

  XEvent report;

  /* Connect to Xserver */
  if ( (*display = XOpenDisplay (disp_name)) == NULL ) {
    fprintf( stderr, "Cannot connect to X server %s\n", XDisplayName(disp_name) );
    return EXIT_FAILURE;
  }

  /* Initialize for graphical display  */
  screen = DefaultScreen( *display );
  disp_width  = DisplayWidth ( *display, screen );
  disp_height = DisplayHeight( *display, screen );
  *win = XCreateSimpleWindow( *display, RootWindow (*display, screen), x, y, width, height,
                border_width, BlackPixel (*display, screen), WhitePixel (*display, screen) );
  XStoreName( *display, *win, window_name );
  *gc = XCreateGC( *display, *win, valuemask, &values ); /* graphics context */
  white = WhitePixel( *display, screen );                /* color value for white */
  black = BlackPixel( *display, screen );                /* color value for black */
  XSetBackground( *display, *gc, white );
  XSetForeground( *display, *gc, black );
  XMapWindow( *display, *win );
  XSync( *display, False );

  /* Get min and max for range of color values
   *   -- assumed to be defined by "white", "black" */
  *min_color = (white > black) ? black : white;
  *max_color = (white > black) ? white : black;

  /* Wait for keyboard input before starting program */
  fprintf( stderr, "Press any key (with focus in display) to start the program\n" );
  fflush( stderr );

  /*  Choose which events we want to handle   */
  XSelectInput( *display, *win, KeyPressMask );

  /* Wait for event */
  XNextEvent( *display, &report );

  return EXIT_SUCCESS;
}



/* Wait for user response before ending program.
 *  -- Also allows user to discover coordinates of points. */
void interact( Display *display, Window *win, uint width, uint height,
               double real_min, double real_max, double imag_min, double imag_max )
{
  double scale_real, scale_imag; 
  XEvent report;
  Window root_return, child_return;
  int root_x_return, root_y_return;
  int win_x_return, win_y_return;
  uint mask_return;

  fprintf( stderr, "\n\n" );
  fprintf( stderr, "Click on a point in the display to get its coordinates\n" );
  fprintf( stderr, "Press any key (with focus in display) to end the program\n" );
  fflush( stderr );

  /* Choose which events we want to handle */
  XSelectInput( display, *win, KeyPressMask | ButtonPressMask );

  /* Compute scaling factors (for processing mouse clicks) */
  scale_real = (double) (real_max - real_min) / (double) width;
  scale_imag = (double) (imag_max - imag_min) / (double) height;

  /* Event loop */
  while (True) {
    XNextEvent( display, &report );

    switch ( report.type ) {

      case ButtonPress:

        XQueryPointer( display, *win, &root_return, &child_return, &root_x_return,
                       &root_y_return, &win_x_return, &win_y_return, &mask_return );
        fprintf( stderr, "coordinates = (%g, %g)\n",
                 real_min + ((double) win_x_return * scale_real),
                 imag_min + ((double) (height-1-win_y_return) * scale_imag));
                   /* height-1-row so y axis displays with larger values at top */
        fflush( stderr );
        break;

      case KeyPress:

        return;
    }
  }
}


/* Function declarations */
int serial_pgm( uint width, uint height, double real_min, double real_max,
                double imag_min, double imag_max, uint maxiter);



/* Draw a row of pixels */
void drawLine( int row, int width, uint *ilines, XDATA *windata )
{
  int col;

  for ( col=0; col<width; ++col ) {
    XSetForeground( windata->display, windata->gc, ilines[width*row+col] );
    XDrawPoint( windata->display, windata->win, windata->gc, col, row );
  }
}



/* Main program */
int main( int argc, char **argv )
{
  int returnval;
  uint maxiter;
  double real_min = -N;
  double real_max =  N;
  double imag_min = -N;
  double imag_max =  N;
  uint width  = NPIXELS;         /* dimensions of display window */
  uint height = NPIXELS;
  double x0 = 0.0, y0 = 0.0;
  double size = 2.0;

  /* Process command-line arguments */
  maxiter = 10000;
  real_min = x0 - size;
  real_max = x0 + size;
  imag_min = y0 - size;
  imag_max = y0 + size;

  MPI_Init( &argc, &argv );
  /* Call serial code */
  returnval = serial_pgm( width, height, real_min, real_max, imag_min, imag_max, maxiter);


  return returnval;
}



/* Program for serial process.
 *  -- Returns EXIT_SUCCESS or EXIT_FAILURE as appropriate. */
int serial_pgm( uint width, uint height, double real_min, double real_max,
                double imag_min, double imag_max, uint maxiter ) {

  XDATA windata;

  ulong min_color, max_color;
  uint col, row;
  int setup_return;

  ulong color;
  double scale_real, scale_imag, scale_color;
  uint count;
  double lengthsq, temp;
  complex z, c;
  uint *ilines;
  double start, finish;

  
  int MPI_Barrier(MPI_Comm comm);
  start = MPI_Wtime();
  
  int processor_id;
  int num_of_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &processor_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_of_proc);
  int initial_position; 
  int final_position;

if (processor_id == 0) {
   
    setup_return = setup(width, height, &windata.display, &windata.win,
                         &windata.gc, &min_color, &max_color);
    if (setup_return != EXIT_SUCCESS) {
        fprintf(stderr, "Unable to initialize display, continuing\n");
    }  

    MPI_Bcast(&windata.win, 1, MPI_UNSIGNED_LONG, 0,MPI_COMM_WORLD);
    MPI_Bcast(&min_color, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_color, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);


    initial_position = 0;
    final_position =(height / num_of_proc);
} else {

    MPI_Bcast(&windata.win, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&min_color, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_color, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    windata.display = XOpenDisplay(NULL);


    windata.gc = XCreateGC(windata.display, windata.win, 0, NULL);
    if(processor_id==num_of_proc-1){
       initial_position = (height / num_of_proc)*processor_id;
       final_position = height;

    }
    else{
    	initial_position = (height / num_of_proc)*processor_id;
	final_position = (height / num_of_proc)*(processor_id+1);
    }

}

  /* Compute factors to scale computational region to window */
  scale_real = (double) (real_max - real_min) / (double) width;
  scale_imag = (double) (imag_max - imag_min) / (double) height; 

  /* Compute factor for color scaling */
  scale_color = (double) (max_color - min_color) / (double) (maxiter - 1);

  /* Allocate memory */
  ilines = calloc( width*height, sizeof(uint) );


  /* Calculate points and draw them a row at a time. */
  for ( row=initial_position; row<final_position; ++row ) {

    /* Scale vertical display coordinates to actual region */
    c.imag = imag_min + ((double) (height-1-row) * scale_imag);
      /* height-1-row so y axis displays with larger values at top */

    for ( col=0; col<width; ++col ) {
      /* Scale horizontal display coordinates to actual region */
      c.real = real_min + ((double) col * scale_real);

      /* Calculate z0,z1,... until divergence or maximum iterations */
      z.real = z.imag = 0.0;
      count = 0;
      do {
        temp   = z.real*z.real - z.imag*z.imag + c.real;
        z.imag = 2.0*z.real*z.imag + c.imag;
        z.real = temp;

        lengthsq = z.real*z.real + z.imag*z.imag;
        ++count;
      } while ( lengthsq<(N*N) && count<maxiter );


      /* Scale color and store */
      color = (ulong) ((count-1) * scale_color) + min_color;

      ilines[width*row+col] = color;
    }


    /* Draw points for this row */
    drawLine( row, width, ilines, &windata );

  }

  /* Be sure everything is written out. */
    XFlush(windata.display);
    MPI_Barrier(MPI_COMM_WORLD);


  int MPI_Barrier(MPI_Comm comm);
  finish = MPI_Wtime();
  printf("Elapsed time = %e seconds\n", finish - start);
  /* Wait for user response, then exit program  */
  if (setup_return==EXIT_SUCCESS) {
    interact( windata.display, &windata.win, width, height,
              real_min, real_max, imag_min, imag_max );
  }

  free(ilines);
  MPI_Finalize();

  return EXIT_SUCCESS;
}
