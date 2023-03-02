import Header from './Header';
import SchemaPage from "../pages/SchemaPage";
// Import and define ThemeProvider at the root of the app
import { createTheme, ThemeProvider } from '@mui/material/styles';


// Global style overrides have changed in MUI5
// Example: https://mui.com/material-ui/customization/theme-components/#global-style-overrides
// Change Table title color
const theme = createTheme({
	components: {
		MUIDataTableToolbar: {
			styleOverrides: {
				titleText: {
					color: "#bb7944",
				},
			},
		},
	}
})


const App = () => {
	return (
		<>
		  <div style={{ marginLeft: "5%", marginRight: "5%", marginTop: "2%", marginBottom: "2%" }}>
				<Header></Header>
				<ThemeProvider theme={theme}>
					<SchemaPage>
					</SchemaPage>
				</ThemeProvider>
		  </div>
		</>
	)
};

export default App;
