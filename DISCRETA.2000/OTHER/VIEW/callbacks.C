// callbacks.C

#include "view.h"




void glxarea_init(Widget w, XtPointer data, XtPointer callData)
{
#ifdef VERBOSE
	printf("in glxarea_init()\n"); fflush(stdout);
#endif
	XtVaGetValues(glx_data->glxarea, 
		XtNwidth, &glx_data->viewWidth, 
		XtNheight, &glx_data->viewHeight, 
		NULL);
#ifdef VERBOSE
	printf("viewWidth=%ld viewHeight=%ld\n", 
		(long int) glx_data->viewWidth, 
		(long int) glx_data->viewHeight);
#endif
	glx_data->glxwin = XtWindow(w);
	if (!glx_data->GLXContext_initialised) {
		if (share_context == NULL) {
			glx_data->cx = glXCreateContext(display, vi, None, True);
	  		share_context = &glx_data->cx;
			}
		else {
			glx_data->cx = glXCreateContext(display, vi, *share_context, True);
			}
		if (glx_data->cx == NULL)
			XtAppError(app_context, "could not create rendering context");
		glx_data->GLXContext_initialised = TRUE;
		}
	glXMakeCurrent(XtDisplay(w), XtWindow(w), glx_data->cx);

	render_init(display_lists_initialised);
	// trackball(curquat, 0.0, 0.0, 0.0, 0.0);
	curquat[0] = appearance.quat[0];
	curquat[1] = appearance.quat[1];
	curquat[2] = appearance.quat[2];
	curquat[3] = appearance.quat[3];
	render_reshape(glx_data, glx_data->viewWidth, glx_data->viewHeight);

	// MGLG->contexts_init_anz++;
  
	// if((MGLG->mit_filearg)&&(! MGLG->start_file_loaded)&&(i2==0)) postLoadFile();

}


void glxarea_resize(Widget w, XtPointer data, XtPointer callData)
{
	GLwDrawingAreaCallbackStruct *resize =
		(GLwDrawingAreaCallbackStruct*) callData;
#ifdef VERBOSE
	printf("in glxarea_resize()\n"); fflush(stdout);
#endif
	int width = resize->width;
	int height = resize->height;
#ifdef VERBOSE
	printf("Width = %d Height = %d\n", width, height); fflush(stdout);
#endif
	render_reshape(glx_data, width, height);
}


void glxarea_draw(Widget w, XtPointer data, XtPointer callData)
{
	// printf("in glxarea_draw()\n"); fflush(stdout);
  postRedisplay();
}

void postRedisplay()
{
	if (!redisplayPending) {
		redisplayID = XtAppAddWorkProc(app_context, handleRedisplay, (void *)0);
		redisplayPending = TRUE;
		}
}

Boolean handleRedisplay(XtPointer callData)
{
	// int i2 = (int) callData;
	if(!redisplayPending)
		return(True);
	renderScene();
#if 0
	if( !zeichen_windows[i2].spinning && !zeichen_windows[i2].pendingAutoHiRes) 
    {
      zeichen_windows[i2].pendingAutoHiRes = 1;
      zeichen_windows[i2].hiResID = XtAppAddTimeOut((*MGLG->app), 100, hiresTimeout, (void*)i2);
    }
#endif
	redisplayPending = FALSE;
	// MGLG->is_drawn[i2] = 1;
	return True;
}




#if 0








void
stopAutoHiRes(int i2)
{
  if(zeichen_windows[i2].pendingAutoHiRes)
    {
      XtRemoveTimeOut(zeichen_windows[i2].hiResID);
      zeichen_windows[i2].pendingAutoHiRes = 0;
    }
}

void
makeHiRes(int i2)
{
  zeichen_windows[i2].moldat.sphereVersion = HI_RES_SPHERE;
  postRedisplay(i2);
}



void
hiresTimeout(XtPointer closure, XtIntervalId *id)
{
  int i2 = (int)closure;
  makeHiRes(i2);
}

void
makeLoRes(int i2)
{
  stopAutoHiRes(i2);
  zeichen_windows[i2].moldat.sphereVersion = LO_RES_SPHERE;
}













#endif


