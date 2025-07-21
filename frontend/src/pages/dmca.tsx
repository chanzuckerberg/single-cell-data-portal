import Head from "next/head";
import Image from "next/image";
import rawCellxgeneLogo from "src/components/common/staticPages/cellxgene.png";
import {
  CommonStyle,
  Layout,
  PrivacyStyle,
} from "src/components/common/staticPages/style";

const DMCA = (): JSX.Element => {
  return (
    <Layout>
      <CommonStyle>
        <PrivacyStyle>
          <Head>
            <title>Digital Millennium Copyright Act (DMCA) Policy</title>
          </Head>
          <header>
            <Image
              data-testid="cellxgene-logo"
              src={rawCellxgeneLogo}
              alt="CELLxGENE logo"
              width="119"
              height="35"
            />
          </header>

          <br />

          <div>
            <h1>Digital Millennium Copyright Act (DMCA) Policy</h1>
            <p>Last updated: February 15, 2025.</p>
            {/* Introduction */}
            <>
              <p>
                We respect the intellectual property rights of others. We
                require that content posted on the Services does not violate the
                intellectual property rights of third parties, including
                copyright and trademark.
              </p>
              <ul className="text-list">
                <li>
                  If you believe that your intellectual property rights have
                  been violated because of content you find on the Services, you
                  can request deletion of your intellectual property rights by
                  following the process described below.
                </li>
                <li>
                  If you have posted content to the Services and someone claims
                  that they own that intellectual property, you will be notified
                  under the process described below.
                </li>
              </ul>
            </>
            {/* Copyright */}
            <>
              <h2>Copyright</h2>
              <p>
                Copyright is a legal right that seeks to protect original works
                of authorship (e.g., books, music, film, art, etc.). Generally,
                copyright protects original expression such as words or images.
                It does not protect facts and ideas, although it may protect the
                original words or images used to describe an idea. Copyright
                also doesn’t protect things like names, titles and slogans;
                however, another legal right called a trademark might protect
                those.
              </p>
              <p>
                For more information on copyright law, you can visit the website
                of the{" "}
                <a href="https://www.copyright.gov/" target="_blank">
                  U.S. Copyright Office
                </a>
                or the
                <a
                  href="http://www.wipo.int/portal/en/index.html"
                  about="_blank"
                >
                  World Intellectual Property Organization (WIPO)
                </a>
                . We can’t provide you with legal advice. Accordingly, please
                consult with an attorney if you have questions about
                intellectual property rights.
              </p>
              <ol>
                <li>
                  <h4>Notification of Infringement.</h4>
                  <p>
                    Chan Zuckerberg Initiative (“CZI”) will respond to notices
                    of alleged copyright infringement that comply with the
                    Digital Millennium Copyright Act (the “DMCA”). In addition,
                    we will promptly terminate without notice the accounts of
                    those determined by us to be “repeat infringers.” If you are
                    a copyright owner or an agent thereof, and you believe that
                    any content hosted on our website or the Services infringes
                    your copyrights, then you may submit a notification pursuant
                    to the DMCA by providing our designated copyright agent
                    (“Designated Agent”) with the following information in
                    writing (please consult your legal counsel or see 17 U.S.C.
                    Section 512(c)(3) to confirm these requirements and your
                    compliance therewith):
                  </p>
                  <ol>
                    <li>
                      A physical or electronic signature of a person authorized
                      to act on behalf of the owner of an exclusive right that
                      is allegedly infringed.
                    </li>
                    <li>
                      Identification of the copyrighted work claimed to have
                      been infringed, or if multiple copyrighted works are
                      covered by a single notification, a representative list of
                      such works at that website.
                    </li>
                    <li>
                      Identification of the material that is claimed to be
                      infringing or to be the subject of infringing activity and
                      that is to be removed or access to which is to be
                      disabled, and information reasonably sufficient to permit
                      us to locate the material. Providing URLs in the body of
                      an email is the best way to help us locate content
                      quickly.
                    </li>
                    <li>
                      Information reasonably sufficient to permit us to contact
                      the complaining party, such as an address, telephone
                      number, and, if available, an electronic mail address at
                      which the complaining party may be contacted.
                    </li>
                    <li>
                      A statement that the complaining party has a good faith
                      belief that use of the material in the manner complained
                      of is not authorized by the copyright owner, its agent, or
                      the law.
                    </li>
                    <li>
                      A statement that the information in the notification is
                      accurate, and under penalty of perjury, that the
                      complaining party is authorized to act on behalf of the
                      owner of an exclusive right that is allegedly infringed.
                    </li>
                  </ol>
                  <p>
                    <strong>
                      Please note that under the DMCA, any person who knowingly
                      materially misrepresents that material or activity is
                      infringing may be subject to liability.
                    </strong>
                  </p>
                  <p>
                    You may submit your notification of alleged copyright
                    infringement by sending it to our “Designated Agent” by mail
                    or e-mail as set forth below. A “Designated Agent” is an
                    individual who has agreed to receive copyright infringement
                    takedown requests for an online service under the DMCA.
                  </p>
                  <p>
                    Please note that we will send a copy of such notices (which
                    will include the personal information you supply in your
                    notice) to the individual who uploaded the allegedly
                    infringing content.
                  </p>
                </li>
                <li>
                  <h4>Counter Notification.</h4>
                  <p>
                    If you elect to send us a counter-notice, after having been
                    notified of a copyright claim submitted to us in accordance
                    with Section 1 above, to be effective it must be a written
                    communication that includes the following (please consult
                    your legal counsel or see 17 U.S.C. Section 512(g)(3) to
                    confirm these requirements):
                  </p>
                  <ol>
                    <li>
                      A physical or electronic signature of the user whose
                      content was removed as a result of the copyright claim.
                    </li>
                    <li>
                      Identification of the material that has been removed or to
                      which access has been disabled and the location at which
                      the material appeared before it was removed or access to
                      it was disabled.
                    </li>
                    <li>
                      A statement under penalty of perjury that the user has a
                      good faith belief that the material was removed or
                      disabled as a result of mistake or misidentification of
                      the material to be removed or disabled.
                    </li>
                    <li>
                      The user’s name, address, and telephone number, and a
                      statement that the user consents to the jurisdiction of
                      Federal District Court for the judicial district in which
                      the address is located, or if the user’s address is
                      outside of the United States, for any judicial district in
                      which we may be found, and that the user will accept
                      service of process from the person who provided
                      notification under Section 1 above or an agent of such
                      person.
                    </li>
                  </ol>
                  <p>
                    <strong>
                      Please note that under the DMCA, any person who knowingly
                      materially misrepresents that material or activity is
                      infringing may be subject to liability.
                    </strong>
                  </p>
                  <p>
                    We only accept counter-notices that meet the requirements
                    set forth above and are received from the email address
                    associated with the account on the Services you used to
                    upload the content within seven (7) business days of our
                    forwarding you the DMCA notice. You may submit your
                    counter-notification by sending it to our Designated Agent
                    by mail or e-mail as set forth below.
                  </p>
                </li>
                <li>
                  <h4>Designated Agent.</h4>
                  <p>
                    Our Designated Agent for notice of claims of copyright
                    infringement can be reached as follows:
                  </p>
                  <br />
                  <br />
                  <address>
                    Attention: Designated Agent Chan Zuckerberg Initiative{" "}
                    <br />
                    2682 Middlefield Road, Suite i <br />
                    Redwood City, CA 94063
                    <br />
                    Attn: General Counsel <br />
                    Email: courtesy copy:{" "}
                    <a href="mailto:privacy@chanzuckerberg.com">
                      privacy@chanzuckerberg.com
                    </a>
                  </address>
                  <p>
                    You acknowledge that if you fail to comply with all of the
                    applicable requirements of this Policy, your DMCA take-down
                    notice may not be valid.
                  </p>
                </li>
              </ol>
            </>
          </div>
        </PrivacyStyle>
      </CommonStyle>
    </Layout>
  );
};

export default DMCA;
