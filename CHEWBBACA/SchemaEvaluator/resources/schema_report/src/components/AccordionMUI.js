import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
// Material-UI ExpansionPanel components
import Accordion from '@mui/material/Accordion';
import AccordionSummary from '@mui/material/AccordionSummary';
import AccordionDetails from '@mui/material/AccordionDetails';
import Typography from '@mui/material/Typography'; 
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';


const AccordionMUI = ({ summaryText, detailsText, expanded }) => {
	return ( 
		<Accordion defaultExpanded={expanded}>
			<AccordionSummary
				expandIcon={<ExpandMoreIcon />}
				aria-controls="panella-content"
				id="panella-header"
			>
				<Typography>{summaryText}</Typography>
			</AccordionSummary>
			<AccordionDetails >
				<div
					//className={classes.mainPaper}
					style={{ width: "100%", height: "100%" }}
				>
					<ReactMarkdown remarkPlugins={[remarkGfm]} children={detailsText} />
				</div>
			</AccordionDetails>
		</Accordion>
	)
};

export default AccordionMUI;
