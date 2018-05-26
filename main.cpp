#include <afxwin.h>  
#include "resource.h"
#define _CRT_SECURE_NO_WARNINGS

#include "Phantom.h" 

class Phantom : public CDialog
{
public:
	Phantom(CWnd* pParent = NULL) : CDialog(Phantom::IDD, pParent)
	{}
	// Dialog Data, name of dialog form
	enum { IDD = IDD_DIALOG1};
protected:
	virtual void DoDataExchange(CDataExchange* pDX) 
	
	{ CDialog::DoDataExchange(pDX); }
	//Called right after constructor. Initialize things here.

	virtual BOOL OnInitDialog()
	{
		CDialog::OnInitDialog();
	
			return true;
	}
public:


	DECLARE_MESSAGE_MAP()

	afx_msg void OnBnClickedinitiliaze();
	afx_msg void OnBnClickedstart();
	afx_msg void OnBnClickedstop();
	afx_msg void OnBnClickedconon();
	afx_msg void OnBnClickedPublisher();
	afx_msg void OnBnClickedCalibarate();
	
	afx_msg void OnBnClickedstoppub();
	afx_msg void OnBnClickedButton1();
	afx_msg void OnBnClickedcamera1();
	afx_msg void OnBnClickedcamera2();
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedCancel();
};
//---------------------------------------------------------------------
class Phantom_winapp : public CWinApp
{
public:
	Phantom_winapp() {  }
public:
	virtual BOOL InitInstance()
	{
		CWinApp::InitInstance();
		Phantom dlg;
		m_pMainWnd = &dlg;
		INT_PTR nResponse = dlg.DoModal();
		return FALSE;
	} //close function

};
//-----------------------For Message Handlers, identifiers and macros---------------------------
//Need a Message Map Macro for both CDialog and CWinApp
BEGIN_MESSAGE_MAP(Phantom, CDialog)   // phantom is owner class name, CDialog is a base class name

	
	ON_BN_CLICKED(CB_initiliaze, &Phantom::OnBnClickedinitiliaze)
	ON_BN_CLICKED(CB_start, &Phantom::OnBnClickedstart)
	ON_BN_CLICKED(CB_stop, &Phantom::OnBnClickedstop)
	ON_BN_CLICKED(CB_con_on, &Phantom::OnBnClickedconon)
	ON_BN_CLICKED(CB_Publisher, &Phantom::OnBnClickedPublisher)
	ON_BN_CLICKED(CB_Calibarate, &Phantom::OnBnClickedCalibarate)
	
	ON_BN_CLICKED(CB_stoppub, &Phantom::OnBnClickedstoppub)
	ON_BN_CLICKED(IDC_BUTTON1, &Phantom::OnBnClickedButton1)
	ON_BN_CLICKED(IDC_camera1, &Phantom::OnBnClickedcamera1)
	ON_BN_CLICKED(IDC_camera2, &Phantom::OnBnClickedcamera2)
	ON_BN_CLICKED(IDOK, &Phantom::OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, &Phantom::OnBnClickedCancel)
END_MESSAGE_MAP()
//-----------------------------------------------------------------------------------------


Phantom_winapp  theApp;  //Starts the Application


void Phantom::OnBnClickedinitiliaze()
{
	// TODO: Add your control notification handler code here
	//if (!(ep = engOpen("")))
	//cout<<"Can't start MATLAB engine"<<endl;

	initiliaze();
	rotation();

	/*
	 * Call engOpen with a NULL string. This starts a MATLAB process 
     * on the current host using the command "matlab".
	 */

	
}

void Phantom::OnBnClickedstart()
{
	// TODO: Add your control notification handler code here

	start();
}


void Phantom::OnBnClickedstop()
{
	// TODO: Add your control notification handler code here

		close();
}


void Phantom::OnBnClickedconon()
{
	// TODO: Add your control notification handler code here
	AllocConsole();
	freopen("CONOUT$", "w", stdout);  //reopen stream with different file or mode 
}



void Phantom::OnBnClickedPublisher()
{
	// TODO: Add your control notification handler code here
	demonstrations();
}


void Phantom::OnBnClickedCalibarate()
{

	// TODO: Add your control notification handler code here
	calibarate();
}



void Phantom::OnBnClickedstoppub()
{
	// TODO: Add your control notification handler code here
	stopDemonstrations();

	/*engEvalString(ep, "close;");*/
}


void Phantom::OnBnClickedButton1()
{
	// TODO: Add your control notification handler code here
}


void Phantom::OnBnClickedcamera1()
{
	


	
}


void Phantom::OnBnClickedcamera2()
{
	// TODO: Add your control notification handler code here
	void camera2();
}


void Phantom::OnBnClickedOk()
{
	// TODO: Add your control notification handler code here
	CDialog::OnOK();
}


void Phantom::OnBnClickedCancel()
{
	// TODO: Add your control notification handler code here
	CDialog::OnCancel();
}
