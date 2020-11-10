import React from "react";
import Header from "./Header";
import SchemaEvaluator from "../components/SchemaEvaluator";

export default function App() {
  return (
    <>
      <div style={{ marginLeft: "5%", marginRight: "5%" }}>
        <Header />
        <SchemaEvaluator />
      </div>
    </>
  );
}
