/**
 * Layout component that queries for data
 * with Gatsby's useStaticQuery component
 *
 * See: https://www.gatsbyjs.org/docs/use-static-query/
 */

import React from "react"
import { useStaticQuery, graphql } from "gatsby"
import { ThemeProvider } from "theme-ui"
import theme from "./theme"

import Header from "./header"
import Footer from "./footer"
import "./layout.css"

const Layout = ({ children }) => {
  const data = useStaticQuery(graphql`
    query SiteTitleQuery {
      site {
        siteMetadata {
          title
        }
      }
    }
  `)

  return (
    <ThemeProvider theme={theme}>
      <Header siteTitle={data.site.siteMetadata.title} />
      <div
        style={{
          margin: `0 auto`,
          maxWidth: 1400,
          minHeight:
            // "100vh", /* for now. better: height minus footer, or go to sticky */
            `calc(100vh - 85px - 60px - 48px)`, // height minus footer, minus header, minus margin between header and "Explore Data"; only way I coudl get footer to level with bottom of screen
          padding: `0 32px`, // since this file seems not to recognize the theme arrays
        }}
      >
        <main>{children}</main>
      </div>
      <Footer />
    </ThemeProvider>
  )
}

export default Layout
