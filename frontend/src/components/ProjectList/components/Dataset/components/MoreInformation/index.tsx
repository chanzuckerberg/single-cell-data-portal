import { FC } from "react";
import { Project } from "src/common/entities";

interface Props {
  links: Project["links"];
}

const MoreInformation: FC<Props> = () => {
  return null;

  // (thuang): Temporarily comment this out until BE has `link.name`
  // return (
  //   <Wrapper>
  //     {links.map(link => {
  //       return (
  //         <StyledAnchor
  //           key={link.url}
  //           href={link.url}
  //           target="_blank"
  //           rel="noopener"
  //         >
  //           {link.name}
  //           <br />
  //         </StyledAnchor>
  //       );
  //     })}
  //   </Wrapper>
  // );
};

export default MoreInformation;
