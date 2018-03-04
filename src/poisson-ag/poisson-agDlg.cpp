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
    , m_bPointsVisible(TRUE)
    , m_bTriangulationVisible(FALSE)
    , m_bDirichletCellsVisible(FALSE)
    , m_bIsolinesVisible(TRUE)
    , m_bFieldLinesVisible(TRUE)
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CPoissonAGDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_PLOT, my_plot);
    DDX_Check(pDX, IDC_CHECK3, m_bPointsVisible);
    DDX_Check(pDX, IDC_CHECK2, m_bTriangulationVisible);
    DDX_Check(pDX, IDC_CHECK1, m_bDirichletCellsVisible);
    DDX_Check(pDX, IDC_CHECK4, m_bIsolinesVisible);
    DDX_Check(pDX, IDC_CHECK5, m_bFieldLinesVisible);
    DDX_Text(pDX, IDC_EDIT1,  data.params->m1_scale);
    DDX_Text(pDX, IDC_EDIT2,  data.params->m1_origin.x);
    DDX_Text(pDX, IDC_EDIT3,  data.params->m1_origin.y);
    DDX_Text(pDX, IDC_EDIT4,  data.params->m1_theta);
    DDX_Text(pDX, IDC_EDIT9,  data.params->m1_qs);
    DDX_Text(pDX, IDC_EDIT10, data.params->m1_qn);
    DDX_Text(pDX, IDC_EDIT5,  data.params->m2_scale);
    DDX_Text(pDX, IDC_EDIT6,  data.params->m2_origin.x);
    DDX_Text(pDX, IDC_EDIT7,  data.params->m2_origin.y);
    DDX_Text(pDX, IDC_EDIT8,  data.params->m2_theta);
    DDX_Text(pDX, IDC_EDIT11, data.params->m2_qs);
    DDX_Text(pDX, IDC_EDIT12, data.params->m2_qn);
    DDX_Text(pDX, IDC_EDIT13, data.params->dx);
    DDX_Text(pDX, IDC_EDIT14, data.params->dy);
}

BEGIN_MESSAGE_MAP(CPoissonAGDlg, CSimulationDialog)
    ON_WM_PAINT()
    ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_BUTTON1, &CPoissonAGDlg::OnBnClickedButton1)
    ON_BN_CLICKED(IDC_BUTTON2, &CPoissonAGDlg::OnBnClickedButton2)
    ON_BN_CLICKED(IDC_CHECK3, &CPoissonAGDlg::OnBnClickedCheck3)
    ON_BN_CLICKED(IDC_CHECK2, &CPoissonAGDlg::OnBnClickedCheck3)
    ON_BN_CLICKED(IDC_CHECK1, &CPoissonAGDlg::OnBnClickedCheck3)
    ON_BN_CLICKED(IDC_CHECK4, &CPoissonAGDlg::OnBnClickedCheck3)
    ON_BN_CLICKED(IDC_CHECK5, &CPoissonAGDlg::OnBnClickedCheck3)
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

    OnBnClickedCheck3();

    my_plot.plot_layer.with(model::make_root_drawable(data.config, {
        data.system_data.dirichlet_cell_plot,
        data.system_data.triangulation_plot,
        data.system_data.system_plot,
        data.isoline_data.plot,
        data.system_data.point_plot
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


void CPoissonAGDlg::OnSimulation()
{
    model::adjust(*data.params, *data.config.world);
    model::update_system_data(*data.params, data.system_data);
    while (m_bWorking)
    {
        Sleep(1000); // stub

        my_plot.RedrawBuffer();
        my_plot.SwapBuffers();

        Invoke([this] () { my_plot.RedrawWindow(); });
    }
    CSimulationDialog::OnSimulation();
}


void CPoissonAGDlg::OnBnClickedButton1()
{
    UpdateData(TRUE);
    StartSimulationThread();
}


void CPoissonAGDlg::OnBnClickedButton2()
{
    UpdateData(TRUE);
    StopSimulationThread();
}


void CPoissonAGDlg::OnBnClickedCheck3()
{
    UpdateData(TRUE);
    data.system_data.point_plot->visible = (m_bPointsVisible == TRUE);
    data.system_data.triangulation_plot->visible = (m_bTriangulationVisible == TRUE);
    data.system_data.dirichlet_cell_plot->visible = (m_bDirichletCellsVisible == TRUE);
    data.isoline_data.plot->visible = (m_bIsolinesVisible == TRUE);
    data.field_line_data.plot->visible = (m_bFieldLinesVisible == TRUE);
}
