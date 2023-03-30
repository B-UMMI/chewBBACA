import MUIDataTable from "mui-datatables";


const DataTable = ({ tableData, tableTitle, tableOptions, tableConditionalFormatting }) => {

	let columnData = (tableData[0].columns).map((column) => {
		return ({ name: column,
				  label: column,
				  options: {
					  filter: true,
					  sort: true,
					  sortThirdClickReset: true,
					  display: true,
				  },
				})
	});

	// add link to locus ID and columns with "_URL" in the name
	if (tableConditionalFormatting) {
		for (let [key, value] of Object.entries(tableConditionalFormatting)) {
			// Add formatting to columns whose name contains the key
			let matchedColumns = (tableData[0].columns).filter(e => e.includes(key));
			for (let columnName of matchedColumns) {
				let columnIndex = tableData[0].columns.indexOf(columnName);
				columnData[columnIndex].options = {...columnData[columnIndex].options, ...value}
			}
		}
	};

	const rowData = tableData[1].rows

	const total_data_table = (
		<MUIDataTable
			title={tableTitle}
			data={rowData}
			columns={columnData}
			options={tableOptions}
		>
		</MUIDataTable>
	  );

	return (
		<div>
			{total_data_table}
		</div>
	)
};

export default DataTable;
