import React, { FC } from "react";
import { Link } from "gatsby";
import { Flex, Box, Image } from "theme-ui";

const HomepageLink = () => {
  return (
    <Link
      to="/"
      style={{
        color: `black`,
        textDecoration: `none`,
      }}
    >
      <Image
        src="https://chanzuckerberg.com/wp-content/themes/czi/img/logo-minified.svg"
        width="55px"
        height="55px"
      />
    </Link>
  );
};

const Footer: FC = () => {
  return (
    <Box
      sx={{
        background: `black`,
        mt: [4],
      }}
    >
      <Flex
        sx={{
          alignItems: "center",
          justifyContent: "space-between",
          fontWeight: 700,
          margin: `0 auto`,
          maxWidth: 1,
          padding: [3],
          paddingLeft: 5,
        }}
      >
        <Box sx={{ width: 112 }}>
          <HomepageLink />
        </Box>
      </Flex>
    </Box>
  );
};

export default Footer;
