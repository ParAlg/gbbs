const vfile = require('to-vfile')
const unified = require('unified')
const markdown = require('remark-parse')
const math = require('remark-math')
const remark2rehype = require('remark-rehype')
const katex = require('rehype-katex')
const stringify = require('rehype-stringify')

unified()
  .use(markdown)
  .use(math)
  .use(remark2rehype)
  .use(katex)
  .use(stringify)
  .process(vfile.readSync('example.md'), function(err, file) {
    if (err) throw err
    console.log(String(file))
  })
