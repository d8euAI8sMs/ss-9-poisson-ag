// poisson-agDlg.cpp : implementation file
//

#include "stdafx.h"
#include "poisson-ag.h"
#include "poisson-agDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CPoissonAGDlg dialog



CPoissonAGDlg::CPoissonAGDlg(CWnd* pParent /*=NULL*/)
    : CSimulationDialog(CPoissonAGDlg::IDD, pParent)
    , data(model::make_model_data())
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CPoissonAGDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_PLOT, my_plot);
}

BEGIN_MESSAGE_MAP(CPoissonAGDlg, CSimulationDialog)
    ON_WM_PAINT()
    ON_WM_QUERYDRAGICON()
END_MESSAGE_MAP()


// CPoissonAGDlg message handlers

BOOL CPoissonAGDlg::OnInitDialog()
{
    CSimulationDialog::OnInitDialog();

    // Set the icon for this dialog.  The framework does this automatically
    //  when the application's main window is not a dialog
    SetIcon(m_hIcon, TRUE);            // Set big icon
    SetIcon(m_hIcon, FALSE);        // Set small icon

    // TODO: Add extra initialization here

    my_plot.plot_layer.with(model::make_root_drawable(data.config, {
        data.system_data.dirichlet_cell_plot,
        data.system_data.triangulation_plot,
        data.system_data.system_plot
    }));

    my_plot.background = plot::palette::brush(0xffffff);
    my_plot.triple_buffered = true;

    my_plot.RedrawBuffer();
    my_plot.SwapBuffers();
    my_plot.RedrawWindow();

    return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CPoissonAGDlg::OnPaint()
{
    if (IsIconic())
    {
        CPaintDC dc(this); // device context for painting

        SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

        // Center icon in client rectangle
        int cxIcon = GetSystemMetrics(SM_CXICON);
        int cyIcon = GetSystemMetrics(SM_CYICON);
        CRect rect;
        GetClientRect(&rect);
        int x = (rect.Width() - cxIcon + 1) / 2;
        int y = (rect.Height() - cyIcon + 1) / 2;

        // Draw the icon
        dc.DrawIcon(x, y, m_hIcon);
    }
    else
    {
        CSimulationDialog::OnPaint();
    }
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CPoissonAGDlg::OnQueryDragIcon()
{
    return static_cast<HCURSOR>(m_hIcon);
}