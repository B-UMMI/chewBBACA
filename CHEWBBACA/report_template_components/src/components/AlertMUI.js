import Alert from '@mui/material/Alert';
import Typography from '@mui/material/Typography';


const AlertMUI = ({ ind, severity, fontSize, children }) => {
	return (
		<Alert key={ind} severity={severity}>
			<Typography key={ind} sx={{ fontSize: fontSize }}>
				{children}
			</Typography>
		</Alert>
	)
}

export default AlertMUI;
