import React from "react";
import ReactDOM from "react-dom/client";
import App from "./components/App";

// get reference to "root" element
const el = document.getElementById('root');
const root = ReactDOM.createRoot(el);
// render the App component
root.render(<App />);
