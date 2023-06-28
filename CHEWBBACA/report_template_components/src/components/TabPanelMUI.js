import { useState } from 'react';
import Box from '@mui/material/Box';
import Tab from '@mui/material/Tab';
import Tabs from '@mui/material/Tabs';
import Paper from '@mui/material/Paper';


function TabPanel(props) {
	const { children, value, index, ...other } = props;
  
	return (
	  <div
		role="tabpanel"
		hidden={value !== index}
		id={`simple-tabpanel-${index}`}
		aria-labelledby={`simple-tab-${index}`}
		{...other}
	  >
		{value === index && (
		  <Box sx={{ p: 3 }}>
			{children}
		  </Box>
		)}
	  </div>
	);
};


function a11yProps(index) {
	return {
		id: `simple-tab-${index}`,
		'aria-controls': `simple-tabpanel-${index}`,
	};
};


const TabPanelMUI = ({ ContentTitles, ContentData }) => {

	const [panel, setPanel] = useState(0);

	const handleChange = (event, newValue) => {
		// where is this value coming from?
		setPanel(newValue);
	};

	// Create Tab component for each title
	const TabTitles = ContentTitles.map((title, index) => {
		return (
			<Tab key={index} label={title} wrapped {...a11yProps(index)} />
		)
	});

	// create TabPanel component for each input data element
	const TabPanelList = ContentData.map((element, index) => {
		return (
			<TabPanel key={index} value={panel} index={index}>
				{element}
			</TabPanel>
		)
	});

	return (
		<Box sx={{ width: "100%" }}>
			<Paper elevation={3}>
				<Box sx={{ borderBottom: 1, borderColor: 'divider' }}>
					<Tabs 
						value={panel} 
						onChange={handleChange} 
						aria-label="basic tabs example" 
						variant="scrollable"
						scrollButtons={true}
						allowScrollButtonsMobile
					>
						{TabTitles}
					</Tabs>
				</Box>
				{TabPanelList}
			</Paper>
		</Box>
	)
};

export default TabPanelMUI;
