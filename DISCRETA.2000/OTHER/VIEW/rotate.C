// rotate.C

#include "view.h"


void startRotation(Widget w, XEvent * event, 
	String * params, Cardinal * num_params)
{
	// printf("in startRotation()\n"); fflush(stdout);
	int x, y;

	if(! graph_loaded) return;
	x = event->xbutton.x;
	y = event->xbutton.y;
	rotate_beginx = x;
	rotate_beginy = y;
	// printf("rotate_beginx=%d rotate_beginy=%d\n", rotate_beginx, rotate_beginy); fflush(stdout);
	stopSpinning();
}

void rotation(Widget w, XEvent * event, 
	String * params, Cardinal * num_params)
{
	// printf("in rotation()\n"); fflush(stdout);
	int x, y;

	if(! graph_loaded) return;
	x = event->xbutton.x;
	y = event->xbutton.y;
	trackball(lastquat,
		(2.0 * rotate_beginx - winWidth) / winWidth,
		(winHeight - 2.0 * rotate_beginy) / winHeight,
		(2.0 * x - winWidth) / winWidth,
		(winHeight - 2.0 * y) / winHeight
		);
	rotate_beginx = x;
	rotate_beginy = y;
	
	if(! spinning) {
		spinning = 1;
		makeLoRes();
		animateID = XtAppAddTimeOut(app_context, 1, animate, (void*) NULL);
		}
}

void animate(XtPointer closure, XtIntervalId *id)
{
	// printf("in animate()\n"); fflush(stdout);
	if(! spinning) {
		if(! pendingAutoHiRes) {
			pendingAutoHiRes = TRUE;
			hiResID = XtAppAddTimeOut(app_context, 100, hiresTimeout, (void*)NULL);
			}
		return;
		}
	//if(!zeichen_windows[j1].redisplayPending)
	add_quats(lastquat, curquat, curquat);
	postRedisplaySpin();
	animateID = XtAppAddTimeOut(app_context, 10, animate, (void*) NULL);
}

void stopSpinning()
{
	if (spinning) {
		spinning = FALSE;
		XtRemoveTimeOut(animateID);
		}
}

void makeHiRes()
{
	sphereVersion = HI_RES_SPHERE;
	postRedisplay();
}

void makeLoRes()
{
	stopAutoHiRes();
	sphereVersion = LO_RES_SPHERE;
}

void stopAutoHiRes()
{
	if (pendingAutoHiRes) {
		XtRemoveTimeOut(hiResID);
		pendingAutoHiRes = FALSE;
		}
}

void postRedisplaySpin()
{
	if(! redisplayPending) {
		redisplayID = XtAppAddWorkProc(app_context, handleRedisplaySpin, (void*) NULL);
		redisplayPending = TRUE;
		}
}

void hiresTimeout(XtPointer closure, XtIntervalId *id)
{
	// int i2 = (int)closure;
	makeHiRes();
}

Boolean handleRedisplaySpin(XtPointer closure)
{
	if(! spinning) {
		stopSpinning();
		redisplayPending = FALSE;
		if(! pendingAutoHiRes) {
			pendingAutoHiRes = TRUE;
			hiResID = XtAppAddTimeOut(app_context, 100, hiresTimeout, (void*) NULL);
			}
		return True;
		}
	redisplayPending = FALSE;
	renderScene();
	if( ! spinning && ! pendingAutoHiRes)  {
		pendingAutoHiRes = TRUE;
		hiResID = XtAppAddTimeOut(app_context, 100, hiresTimeout, (void*) NULL);
		}
  
	return True;
}



