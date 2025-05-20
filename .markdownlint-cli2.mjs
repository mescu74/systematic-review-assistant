/**
 * Configured to align with Markdown style used by most frontier models and also
 */
import configOptions, { init } from "@github/markdownlint-github";

const overriddenOptions = init({
	"single-h1": false,
	"ul-style": "dash",
	"emphasis-style": "asterisk",
	"required-headings": false,
	"list-marker-space": {
		ul_multi: 1,
		ul_single: 1,
		start_indented: false,
	},
	"ul-indent": {
		indent: 4,
	},
	"ol-indent": {
		indent: 4,
	},
	"line-length": {
		code_blocks: false,
		tables: false,
		headings: true,
		heading_line_length: 100,
		line_length: 800,
	},
	"no-inline-html": false,
	"no-emphasis-as-heading": false,
	"blanks-around-lists": false,
	"blanks-around-headers": false,
	"blanks-around-fences": false,
	"blanks-around-lists": false,
	"blanks-around-headers": false,
	"blanks-around-fences": false,
	//"blanks-around-code-blocks": false,
	"blanks-around-tables": true,
	//"blanks-around-images": false,
	//"blanks-around-links": false,
	"blanks-around-blockquotes": true,
	"reference-links-images": {
		shortcut_syntax: false,
	},
	"no-inline-html": false,
	"no-emphasis-as-heading": false,
});

/** @type {import("markdownlint-cli2").Options} */
const options = {
	config: overriddenOptions,
	fix: true,
	ignores: ["**/node_modules", "**/.*"],
	gitignore: true,
	customRules: ["@github/markdownlint-github"],
	outputFormatters: [
		["markdownlint-cli2-formatter-pretty", { appendLink: true }], // ensures the error message includes a link to the rule documentation
	],
};
export default options;
