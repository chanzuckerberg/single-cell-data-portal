import React from "react"
import { Link } from "gatsby"

import Layout from "../components/layout"
import SEO from "../components/seo"

const SecondPage = () => (
  <Layout>
    <SEO title="Projects" />
    <h1>Foo</h1>
    <p>Bar</p>
    <Link to="/">url to home</Link>
  </Layout>
)

export default SecondPage
