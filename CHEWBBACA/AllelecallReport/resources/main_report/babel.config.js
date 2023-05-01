'use strict'
/**
babel.config.js with useful plugins. 
*/
module.exports = function(api) {
  api.cache(true);

  const presets = [
                    [
                      "@babel/preset-env", {
                        "targets": {
                          "esmodules": true,
                          "node":true
                        }
                      }
                    ],
					// add to enable support for experimental 'jsx'
					["@babel/preset-react", {"runtime": "automatic"}],
                  ];
  const plugins = [
	  ["@babel/plugin-proposal-class-properties"]
  ];

  return {
    presets,
    plugins
  }
}
