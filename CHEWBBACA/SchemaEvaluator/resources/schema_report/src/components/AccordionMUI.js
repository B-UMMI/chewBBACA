// Material-UI ExpansionPanel components
import Accordion from '@mui/material/Accordion';
import AccordionSummary from '@mui/material/AccordionSummary';
import AccordionDetails from '@mui/material/AccordionDetails';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import Divider from '@mui/material/Divider';

const AccordionMUI = ({ summaryText, detailsText, expanded }) => {
	return ( 
		<Accordion defaultExpanded={expanded}>
			<AccordionSummary
				expandIcon={<ExpandMoreIcon />}
				aria-controls="panella-content"
				id="panella-header"
			>
				{summaryText}
			</AccordionSummary>
			<Divider></Divider>
			<AccordionDetails >
				<div
					style={{ width: "100%", height: "100%" }}
				>
					{detailsText}
				</div>
			</AccordionDetails>
		</Accordion>
	)
};

export default AccordionMUI;
