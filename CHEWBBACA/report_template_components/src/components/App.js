import ReportPage from "../pages/ReportPage";
// Import and define ThemeProvider at the root of the app
import { createTheme, ThemeProvider } from '@mui/material/styles';


// Global style overrides have changed in MUI5
// Example: https://mui.com/material-ui/customization/theme-components/#global-style-overrides
// Change Table title color
const theme = createTheme({
	components: {
		// Override style in mui-datatables table title
		MUIDataTableToolbar: {
			styleOverrides: {
				titleText: {
					color: "#bb7944",
					//fontWeight: 'bold'
				},
			},
		},
		// Override style in mui-datatables table header cells
		MUIDataTableHeadCell: {
			styleOverrides: {
				data: {
					fontWeight: 'bold',
					textTransform: 'none',
					textAlign: 'left'
				}
			}
		},
		// Override style in mui-datatables table body cells
		// MUIDataTableBodyCell: {
		// 	styleOverrides: {
		// 		root: {
		// 			fontWeight: "bold",
		// 		}
		// 	}
		// },
		// Override style in MUI tabs
		MuiTab: {
			styleOverrides: {
				root: {
					fontWeight: 'bold',
					textTransform: 'none'
				}
			}
		}
	}
})


const App = () => {
	return (
		<div style={{ marginLeft: "5%", marginRight: "5%", marginTop: "2%", marginBottom: "2%" }}>
			<ThemeProvider theme={theme}>
				<ReportPage>
				</ReportPage>
			</ThemeProvider>
		</div>
	)
};

export default App;
