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
    , m_nIsolineCount(100)
    , m_fpIsolineDelta(0.1)
    , m_nFieldLineCount(100)
    , simmode(simmode_default)
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
    this->params = *data.params;
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
    DDX_Text(pDX, IDC_EDIT1, params.w);
    DDX_Text(pDX, IDC_EDIT5, params.h);
    DDX_Text(pDX, IDC_EDIT2, params.r);
    DDX_Text(pDX, IDC_EDIT3, params.dr.y);
    DDX_Text(pDX, IDC_EDIT10, params.dr2.x);
    DDX_Text(pDX, IDC_EDIT11, params.dr2.y);
    DDX_Text(pDX, IDC_EDIT4, params.s);
    DDX_Text(pDX, IDC_EDIT9, params.ds);
    DDX_Text(pDX, IDC_EDIT13, params.q0);
    DDX_Text(pDX, IDC_EDIT14, params.q02);
    DDX_Text(pDX, IDC_EDIT21, params.eps);
    DDX_Text(pDX, IDC_EDIT6, params.b.x);
    DDX_Text(pDX, IDC_EDIT7, params.b.y);
    DDX_Text(pDX, IDC_EDIT15, m_nIsolineCount);
    DDX_Text(pDX, IDC_EDIT16, m_fpIsolineDelta);
    DDX_Text(pDX, IDC_EDIT17, m_nFieldLineCount);
    DDX_Control(pDX, IDC_EDIT12, udata.U1);
    DDX_Control(pDX, IDC_EDIT19, udata.U2);
    DDX_Control(pDX, IDC_EDIT20, udata.Ud);
    DDX_Control(pDX, IDC_EDIT22, udata.d1x);
    DDX_Control(pDX, IDC_EDIT24, udata.d1y);
    DDX_Control(pDX, IDC_EDIT23, udata.d2x);
    DDX_Control(pDX, IDC_EDIT25, udata.d2y);
    DDX_Control(pDX, IDC_COMBO1, m_simmode);
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
    ON_CBN_SELCHANGE(IDC_COMBO1, &CPoissonAGDlg::OnCbnSelchangeCombo1)
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
        data.field_line_data.plot,
        data.system_data.point_plot
    }));

    my_plot.background = plot::palette::brush();
    my_plot.triple_buffered = true;

    my_plot.RedrawBuffer();
    my_plot.SwapBuffers();
    my_plot.RedrawWindow();

    m_simmode.AddString(_T("Default"));
    m_simmode.AddString(_T("Reverse"));
    m_simmode.AddString(_T("d1 calculation"));
    m_simmode.AddString(_T("d2 calculation"));
    m_simmode.SetCurSel(0);

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
    *data.params = params;

    switch (simmode)
    {
    case simmode_rev:
        data.params->b = { -params.b.x, -params.b.y };
        data.params->q0 = params.q02;
        data.params->q02 = params.q0;
        data.params->dr = params.dr2;
        data.params->dr2 = params.dr;
        data.params->has_q0 = TRUE;
        data.params->has_q02 = FALSE;
        break;
    case simmode_dipole2_calc:
        data.params->b = { -params.b.x, -params.b.y };
        data.params->q0 = params.q02;
        data.params->q02 = params.q0;
        data.params->dr = params.dr2;
        data.params->dr2 = params.dr;
        data.params->has_q0 = TRUE;
        data.params->has_q02 = FALSE;
        data.params->has_1 = TRUE;
        data.params->has_2 = FALSE;
        break;
    case simmode_dipole1_calc:
        data.params->has_1 = TRUE;
        data.params->has_2 = FALSE;
        break;
    default:
    case simmode_default:
        break;
    }

    model::adjust(*data.params, *data.config.world);
    model::update_system_data(*data.params, data.system_data);
    model::finel_galerkin g(*data.params, *data.system_data.geometry, data.system_data.mesh);
    while (m_bWorking)
    {
        g.next(*data.system_data.data);

        model::find_isolines
        (
            *data.system_data.mesh,
            *data.system_data.data,
            m_fpIsolineDelta,
            m_nIsolineCount,
            *data.isoline_data.data
        );

        std::vector < geom::mesh::idx_t > hints =
            data.system_data.geometry->hints;

        if (hints.size() > m_nFieldLineCount)
        {
            size_t n = (size_t) std::floor(hints.size() / m_nFieldLineCount + 1);
            if (n > 1)
            {
                hints.resize(hints.size() / n);
                for (size_t i = 0; i < hints.size(); ++i)
                {
                    hints[i] = data.system_data.geometry->hints[i * n];
                }
            }
        }

        model::field_line_finder fnd(
            data.system_data.mesh,
            *data.system_data.data,
            hints);
        fnd.find(*data.field_line_data.data);

        switch (simmode)
        {
        case simmode_default:
        {
            auto i = data.system_data.mesh->find_nearest(data.params->dr2 * data.params->r + data.params->b);
            double u = data.system_data.data->at(i) * data.params->q02;
            CString fmt; fmt.Format(_T("%.9f"), u);
            udata.U1.SetWindowText(fmt);
        }
            break;
        case simmode_rev:
        {
            auto i = data.system_data.mesh->find_nearest(data.params->dr2 * data.params->r + data.params->b);
            double u = data.system_data.data->at(i) * data.params->q02;
            CString fmt; fmt.Format(_T("%.9f"), u);
            udata.U2.SetWindowText(fmt);
        }
            break;
        case simmode_dipole1_calc:
        {
            udata.d1 = DipoleAt(data.params->b, { 0, 0 });
            CString fmt; fmt.Format(_T("%.6f"), udata.d1.x);
            udata.d1x.SetWindowText(fmt);
            fmt.Format(_T("%.6f"), udata.d1.y);
            udata.d1y.SetWindowText(fmt);
            if (udata.d2 != geom::point2d_t(0, 0))
            {
                auto o = data.params->b / data.params->b.length();
                double u = udata.d1.x * udata.d2.x + udata.d1.y * udata.d2.y;
                u += 3 * (udata.d1.x * o.x + udata.d1.y * o.y) * (udata.d2.x * o.x + udata.d2.y * o.y);
                u /= data.params->b.length() * data.params->b.length() * data.params->b.length();
                fmt.Format(_T("%.5e"), u);
                udata.Ud.SetWindowText(fmt);
            }
        }
            break;
        case simmode_dipole2_calc:
        {
            udata.d2 = DipoleAt(data.params->b, { 0, 0 });
            udata.d2 = { -udata.d2.x, -udata.d2.y };
            CString fmt; fmt.Format(_T("%.6f"), udata.d2.x);
            udata.d2x.SetWindowText(fmt);
            fmt.Format(_T("%.6f"), udata.d2.y);
            udata.d2y.SetWindowText(fmt);
            if (udata.d1 != geom::point2d_t(0, 0))
            {
                auto o = data.params->b / data.params->b.length();
                double u = udata.d1.x * udata.d2.x + udata.d1.y * udata.d2.y;
                u += 3 * (udata.d1.x * o.x + udata.d1.y * o.y) * (udata.d2.x * o.x + udata.d2.y * o.y);
                u /= data.params->b.length() * data.params->b.length() * data.params->b.length();
                fmt.Format(_T("%.5e"), u);
                udata.Ud.SetWindowText(fmt);
            }
        }
            break;
        default:
            break;
        }

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


geom::point2d_t CPoissonAGDlg::DipoleAt(geom::point2d_t p0, geom::point2d_t o)
{
    auto v1 = data.system_data.mesh->find_nearest(p0);
    if (v1 == SIZE_T_MAX) return { 0, 0 };

    geom::point2d_t p = { 0, 0 };

    size_t i = 0;

    for each (auto v2 in data.system_data.mesh->vertices()[v1].neighbor_vertices)
    {
        auto p01 = data.system_data.mesh->point_at(v1);
        auto p02 = data.system_data.mesh->point_at(v2);

        auto p1 = (p01 - o) / math::norm(p01 - o);
        auto p2 = (p02 - o) / math::norm(p02 - o);

        double dn = p2.y * p1.x - p2.x * p1.y;

        geom::point2d_t d = {
            (p1.length() * p1.length() * data.system_data.data->at(v1) * p2.y -
            p2.length() * p2.length() * data.system_data.data->at(v2) * p1.y) / dn,
            -(p1.length() * p1.length() * data.system_data.data->at(v1) * p2.x -
            p2.length() * p2.length() * data.system_data.data->at(v2) * p1.x) / dn
        };

        if (!isfinite(d.x) || !isfinite(d.y)) continue;

        p += d;
        ++i;
    }

    p = p / i;

    return p;
}


void CPoissonAGDlg::OnCbnSelchangeCombo1()
{
    simmode = (_simmode) m_simmode.GetCurSel();
}
