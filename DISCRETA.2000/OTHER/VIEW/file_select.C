// file_select.C



#include "view.h"

#define VERBOSE

FILE_SELECT_DIALOG *fsel_open(
	Display *display, char *app_name, int app_argc, char **app_argv, 
	char *filter_text, 
	char *list_text, 
	char *select_text, 
	char *pattern, 
	int f_name_used,
	char *f_name, 
	void (* do_func)(FILE_SELECT_DIALOG *fsel))
{
	FILE_SELECT_DIALOG *fsel = NULL;
	char *p;
	Arg al[64];                    /* Arg List */
	register int ac = 0;           /* Arg Count */
	
	fsel = (FILE_SELECT_DIALOG *) malloc(sizeof(FILE_SELECT_DIALOG));
	if (fsel == NULL) {
		XtAppError(app_context, "no memory for FILE_SELECT_DIALOG");
		}
	strcpy(fsel->filter_text, filter_text);
	strcpy(fsel->list_text, list_text);
	strcpy(fsel->select_text, select_text);
	strcpy(fsel->pattern, pattern);
	fsel->f_name_used = f_name_used;
	strcpy(fsel->f_name, f_name);
	fsel->do_func = do_func;
	if (f_name_used) {
		strcpy(fsel->directory, f_name);
		p = strrchr(fsel->directory, '/');
		if (p) {
			p++;
			*p = 0;
			}
		else
			strcpy(fsel->directory, ".");
		}
	else {
		strcpy(fsel->directory, ".");
		}
#ifdef VERBOSE
	printf("fsel structure initialized\n");
	printf("directory=%s\n", fsel->directory);
	printf("f_name=%s\n", fsel->f_name);
	printf("pattern=%s\n", fsel->pattern);
	fflush(stdout);
#endif
	{
	Arg wargs[32];
	int n = 0;

#if 1
	fsel->shell = XtAppCreateShell ( app_name, "FileselectShell", applicationShellWidgetClass, display, al, ac );
#else
	fsel->shell = XmCreateDialogShell(s1->shell1, "FileselectShell", NULL, 0);
#endif

#if 0
	fsel->shell = XtCreateApplicationShell("FileselectShell", 
		transientShellWidgetClass, /* oder: topLevelShellWidgetClass */
		NULL, 0);
	printf("application shell created\n"); fflush(stdout);
	fsel->shell = XtCreateWidget("FileselectShell", 
		xmDialogShellWidgetClass, 
		NULL, 0);
	printf("dialog shell created\n"); fflush(stdout);
#endif
	
	
#if 1
	if (fsel->f_name_used) {
		XmString text = (XmString)NULL;
		text = XmStringCreate(fsel->f_name, XmSTRING_DEFAULT_CHARSET);
		XtSetArg(wargs[n], XmNdirSpec/*XmNtextString*/, text); n++;
		}



	XmString filter_text = (XmString)NULL;
	XmString list_text = (XmString)NULL;
	XmString select_text = (XmString)NULL;
	filter_text = XmStringCreate(fsel->filter_text, XmSTRING_DEFAULT_CHARSET);
	list_text = XmStringCreate(fsel->list_text, XmSTRING_DEFAULT_CHARSET);
	select_text = XmStringCreate(fsel->select_text, XmSTRING_DEFAULT_CHARSET);
	XtSetArg(wargs[n], XmNfilterLabelString, filter_text); n++;
	XtSetArg(wargs[n], XmNlistLabelString, list_text); n++;
	XtSetArg(wargs[n], XmNselectionLabelString, select_text); n++;



	XmString directory = (XmString)NULL;
	directory = XmStringCreate(fsel->directory, XmSTRING_DEFAULT_CHARSET);
	XtSetArg(wargs[n], XmNdirectory, directory); n++;

#if 1
	XmString pattern = (XmString)NULL;
	pattern = XmStringCreate(fsel->pattern, XmSTRING_DEFAULT_CHARSET);
	XtSetArg(wargs[n], XmNpattern, pattern); n++;
#endif

#if 1
	char dp[10000];
	XmString dp_string = (XmString)NULL;
	strcpy(dp, fsel->directory);
	strcat(dp, "/");
	strcat(dp, fsel->pattern);
	XtSetArg(wargs[n], XmNdirMask, dp_string); n++;
#endif
#endif


	/*XtSetValues(shell, wargs, n);*/
	fsel->fsel_box = XtCreateManagedWidget("fselbox", 
		xmFileSelectionBoxWidgetClass, 
		fsel->shell, wargs, n);
#ifdef VERBOSE
	printf("fsel box created\n"); fflush(stdout);
#endif
	XtUnmanageChild(XmFileSelectionBoxGetChild(fsel->fsel_box, XmDIALOG_HELP_BUTTON));


	/* Callbacks: */
	XtAddCallback(fsel->fsel_box, XmNapplyCallback, cb_file_select, fsel);
	XtAddCallback(fsel->fsel_box, XmNcancelCallback, cb_file_select, fsel);
	XtAddCallback(fsel->fsel_box, XmNnoMatchCallback, cb_file_select, fsel);
	XtAddCallback(fsel->fsel_box, XmNokCallback, cb_file_select, fsel);
	
	fsel->f_realized = FALSE;
	// XtRealizeWidget(fsel->shell);

	}

	return fsel;
}

void fsel_realize(FILE_SELECT_DIALOG *fsel, int mode)
{
	if (fsel->f_realized)
		return;
	XtRealizeWidget(fsel->shell);
	fsel->f_realized = TRUE;
	fsel->mode = mode;
}

void fsel_unrealize(FILE_SELECT_DIALOG *fsel)
{
	if (!fsel->f_realized)
		return;
	XtUnrealizeWidget(fsel->shell);
	fsel->f_realized = FALSE;
}


#if 0
INT x11_fileselect(
	BYTE *filter_text, 
	BYTE *list_text, 
	BYTE *select_text, 
	BYTE *mask, 
	INT *f_name_used,
	BYTE *f_name, 
	void *data,
	void (*do_func)(void *data))
/* 
 * filter_text: Titelzeile, z.B (default bei NIL) 'Datei Filter', 
 * list_text: oberhalb der Files Listbox, z.B. 'PS - Datei' (default Datei), 
 * select_text: oberhalb der Dateiauswahlzeile, z.B. (default) 'Auswahl' 
 * mask: (default) '*.*' 
 */
{
	FILE_SELECT_DIALOG *fsel;

	fsel = (FILE_SELECT_DIALOG *) malloc(sizeof(FILE_SELECT_DIALOG));
	if (fsel == NIL) {
		printf("x11_fileselect() no memory\n");
		return ERROR;
		}
	fsel->data = data;
	fsel->do_func = do_func;
	if (filter_text)
		fsel->filter_text = filter_text;
	else
		fsel->filter_text = "Datei Filter";
	if (list_text)
		fsel->list_text = list_text;
	else
		fsel->list_text = "Datei";
	if (fsel->select_text)
		fsel->select_text = select_text;
	else
		fsel->select_text = "Auswahl:";
	if (mask)
		fsel->mask = mask;
	else
		fsel->mask = "*.*";
	fsel->f_name_used = f_name_used;
	fsel->f_name = f_name;
	{
	Arg wargs[32];
	INT n;
	XmString mask = (XmString)NULL;
	XmString text = (XmString)NULL;
	XmString filter_text = (XmString)NULL;
	XmString list_text = (XmString)NULL;
	XmString select_text = (XmString)NULL;
	
	fsel->shell = XtCreateApplicationShell("FileselectShell", 
		transientShellWidgetClass, /* oder: topLevelShellWidgetClass */
		NULL, 0);
	mask = XmStringCreate(tstr(fsel->mask), 
		XmSTRING_DEFAULT_CHARSET);
	if (*fsel->f_name_used) {
		text = XmStringCreate(tstr(fsel->f_name), 
			XmSTRING_DEFAULT_CHARSET);
		}
	filter_text = XmStringCreate(fsel->filter_text, 
		XmSTRING_DEFAULT_CHARSET);
	list_text = XmStringCreate(fsel->list_text, 
		XmSTRING_DEFAULT_CHARSET);
	select_text = XmStringCreate(fsel->select_text, 
		XmSTRING_DEFAULT_CHARSET);
	n = 0;
	XtSetArg(wargs[n], XmNdirMask, mask); n++;
	XtSetArg(wargs[n], XmNdirSpec/*XmNtextString*/, text); n++;
	XtSetArg(wargs[n], XmNfilterLabelString, filter_text); n++;
	XtSetArg(wargs[n], XmNlistLabelString, list_text); n++;
	XtSetArg(wargs[n], XmNselectionLabelString, select_text); n++;
	/*XtSetValues(shell, wargs, n);*/
	fsel->fsel_box = XtCreateManagedWidget("fselbox", 
		xmFileSelectionBoxWidgetClass, 
		fsel->shell, wargs, n);
	
	/* Callbacks: */
	XtAddCallback(fsel->fsel_box, XmNapplyCallback, cb_file_select, fsel);
	XtAddCallback(fsel->fsel_box, XmNcancelCallback, cb_file_select, fsel);
	XtAddCallback(fsel->fsel_box, XmNnoMatchCallback, cb_file_select, fsel);
	XtAddCallback(fsel->fsel_box, XmNokCallback, cb_file_select, fsel);
	
	XtRealizeWidget(fsel->shell);

	}
	return OK;
}
#endif

void cb_file_select(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	FILE_SELECT_DIALOG *fsel = (FILE_SELECT_DIALOG *) client_data;
	XmSelectionBoxCallbackStruct *call_data = 
		(XmSelectionBoxCallbackStruct *) xt_call_data;
	int n;
	Arg wargs[32];
	XmString text = (XmString)NULL;
	char *s = NULL;
	char str[256];
	
	n=0;
	XtSetArg(wargs[n], XmNdirSpec, &text); n++;
	XtGetValues(fsel->fsel_box, wargs, n);
	XmStringGetLtoR(text, XmSTRING_DEFAULT_CHARSET, &s);
	if (call_data->reason == XmCR_OK) {
		printf("cb_file_select()|reason: XmCR_OK\n");
		sprintf(str, "chosen: %s", s);
		printf("%s\n", str);
		/* strcpy(fsel->fname, s); */
		strcpy(fsel->f_name, s);
		fsel->f_name_used = TRUE;
		if (fsel->do_func == NULL) {
			printf("cb_file_select()|fsel->do_func == NIL\n");
			return;
			}
		(*fsel->do_func)(fsel);
		}
	else if (call_data->reason == XmCR_APPLY) {
		printf("cb_file_select()|reason: XmCR_APPLY\n");
		}
	else if (call_data->reason == XmCR_NO_MATCH) {
		printf("cb_file_select()|reason: XmCR_NO_MATCH\n");
		}
	else if (call_data->reason == XmCR_CANCEL) {
		printf("cb_file_select()|reason: XmCR_CANCEL\n");
		// XtUnrealizeWidget(fsel->shell);
		// free(fsel);
		fsel_unrealize(fsel);
		}
	else if (call_data->reason == XmCR_HELP) {
		printf("cb_file_select()|reason: XmCR_HELP\n");
		}
}


