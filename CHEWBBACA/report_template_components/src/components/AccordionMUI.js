// Material-UI ExpansionPanel components
import Accordion from '@mui/material/Accordion';
import AccordionSummary from '@mui/material/AccordionSummary';
import AccordionDetails from '@mui/material/AccordionDetails';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import Divider from '@mui/material/Divider';
import Stack from '@mui/material/Stack';


const AccordionMUI = ({ summaryText, detailsData, expanded, alerts }) => {
	return (
		<div>
			<Accordion defaultExpanded={expanded}>
				<AccordionSummary
					expandIcon={<ExpandMoreIcon />}
					aria-controls="panella-content"
					id="panella-header"
				>
					{summaryText}
				</AccordionSummary>
				<Divider></Divider>
				<AccordionDetails>
					<Stack sx={{ width: "100%", marginTop: "10px", marginBottom: "20px" }} spacing={1}>
						{alerts}
					</Stack>
					<div style={{ width: "100%", height: "100%" }}>
						{detailsData}
					</div>
				</AccordionDetails>
			</Accordion>
		</div>
	)
};

export default AccordionMUI;
