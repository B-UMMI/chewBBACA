// Material-UI components
import Accordion from '@mui/material/Accordion';
import AccordionSummary from '@mui/material/AccordionSummary';
import AccordionDetails from '@mui/material/AccordionDetails';
import Typography from '@mui/material/Typography'; 
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import Divider from '@mui/material/Divider';


const AccordionMUI = ({ summaryText, detailsData, expanded }) => {
	return (
		<Accordion defaultExpanded={expanded}>
			<AccordionSummary
				expandIcon={<ExpandMoreIcon />}
				aria-controls="panella-content"
				id="panella-header"
			>
				<Typography>{summaryText}</Typography>
			</AccordionSummary>
			<Divider />
			<AccordionDetails >
				<div style={{ width: "100%", height: "100%" }} >
					{detailsData}
				</div>
			</AccordionDetails>
		</Accordion>
	)
};

export default AccordionMUI;
