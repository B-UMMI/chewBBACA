import Alert from '@mui/material/Alert';
import Typography from '@mui/material/Typography';


const AlertMUI = ({ severity, fontSize, children }) => {
	return (
		<Alert severity={severity}>
			<Typography sx={{ fontSize: fontSize }}>
				{children}
			</Typography>
		</Alert>
	)
}

export default AlertMUI;
