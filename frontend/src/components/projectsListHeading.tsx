import React, { FC } from "react";
import { Flex, Box } from "theme-ui";
import { SystemStyleObject } from "@styled-system/css";

const boxBaseStyle: SystemStyleObject = {
  fontSize: [2],
  boxSizing: "border-box",
  flexGrow: 0,
  flexShrink: 0,
  width: `16%`, // Narrower so Project Name has more width; account for right margin of "6" in array = 48px
  paddingRight: 6,
  overflow: "hidden", // Or flex might break
  listStyle: "none",
  fontWeight: 700,
};

const StyledBox: FC = props => {
  return <Box {...props} sx={boxBaseStyle} />;
};

const CellsBox: FC = props => {
  const style: SystemStyleObject = {
    ...boxBaseStyle,
    paddingRight: 7,
    textAlign: "right", // since it holds numbers
  };

  return <Box {...props} sx={style} />;
};

const ProjectBox: FC = props => {
  const style: SystemStyleObject = {
    ...boxBaseStyle,
    width: `36%`, // Narrower so Project Name has more width; account for right margin of "6" in array = 48px
    paddingRight: "unset",
  };

  return <Box {...props} sx={style} />;
};

const ProjectsListHeading: FC = () => {
  return (
    <Flex sx={{ mb: [4] }}>
      <StyledBox>Organ</StyledBox>
      <StyledBox>Assay</StyledBox>
      <StyledBox>Species</StyledBox>
      <CellsBox>Cells</CellsBox>
      <ProjectBox>Project name</ProjectBox>
    </Flex>
  );
};

export default ProjectsListHeading;
